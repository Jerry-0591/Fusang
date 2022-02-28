import os
import re
import sys
import numpy as np
import functools
import argparse
from Bio import AlignIO
from itertools import combinations

import math
import random
import pandas as pd
import tensorflow as tf
from tensorflow.keras import layers, models, Sequential, regularizers

from ete3 import Tree
from pandas.core.frame import DataFrame

global org_seq, comb_of_id, dl_model, dl_predict, len_of_msa, dic_for_leave_node_comb_name, start_end_list
global taxa_num, leave_node_id, leave_node_name, leave_node_comb_id, leave_node_comb_name, internal_node_name_pool

from multiprocessing import Process, Pool
import multiprocessing

#import logging
#logging.basicConfig(level = logging.INFO, filename='new.log', filemode='a')


def comb_math(n,m):
    return math.factorial(n)//(math.factorial(n-m)*math.factorial(m))


def nlargest_indices(arr, n):
    uniques = np.unique(arr)
    threshold = uniques[-n]
    return np.where(arr >= threshold)


def get_quartet_ID(quartet):
    return "".join(sorted(quartet))


def tree_from_quartet(quartet):
    root = Tree()
    root.name = "internal_node_0"
    left = root.add_child(name="internal_node_1")
    left.add_child(name=quartet[0])
    left.add_child(name=quartet[1])
    right = root.add_child(name="internal_node_2")
    right.add_child(name=quartet[2])
    right.add_child(name=quartet[3])
    for desc in root.iter_descendants():
        desc.dist = 0
    return root


def get_topology_ID(quartet):
    return get_quartet_ID(quartet[0:2]) + get_quartet_ID(quartet[2:4])


def get_current_topology_id(quart_key, cluster_1, cluster_2):
    ans = []
    a1 = quart_key.index(cluster_1)
    a2 = quart_key.index(cluster_2)
    ans.append(str(a1))
    ans.append(str(a2))
    ans = set(ans)
    if ans == {'0', '1'} or ans == {'2', '3'}:
        return 0
    elif ans == {'0', '2'} or ans == {'1', '3'}:
        return 1
    elif ans == {'0', '3'} or ans == {'1', '2'}:
        return 2
    else:
        print('Error of function get_current_topology_id, exit the program')
        sys.exit(0)


def judge_tree_score(tree, quart_distribution, new_addition_taxa, dic_for_leave_node_comb_name):
    '''
        parameter 
        tree: a candidate tree, can be any taxas
        quart_distribution: the prob distribution of the topology of every 4-taxa
    '''
    crt_tree = tree.copy("newick")
    leaves = crt_tree.get_leaves()
    
    leaves = [ele.name for ele in leaves]
    total_quarts = list(combinations(leaves, 4))
    quarts = []
    for ele in total_quarts:
        if new_addition_taxa in ele:
            quarts.append(ele)
    
    total_quart_score = 0
    
    for quart in quarts:
        crt_tree = tree.copy("newick")
        try:
            crt_tree.prune(list(quart))
        except:
            print('Error of pruning 4 taxa from current tree, the current tree is:')
            print(crt_tree)
            sys.exit(0)
        
        quart_key = "".join(sorted(list(quart)))
        #quart_topo_id = leave_node_comb_name.index(quart_key)

        quart_topo_id = dic_for_leave_node_comb_name[quart_key]

        quart_topo_distribution = quart_distribution[quart_topo_id]
        
        # judge current tree belongs to which topology
        tmp = re.findall("\([\s\S]\,[\s\S]\)", crt_tree.write(format=9))[0]
        topology_id = get_current_topology_id(quart_key, tmp[1], tmp[3])
        
        total_quart_score += np.log(quart_topo_distribution[topology_id]+1e-200)
    
    return total_quart_score


def get_modify_tree(tmp_tree, edge_0, edge_1, new_add_node_name):
    '''
        add a new leave node between edge_0 and edge_1
        default: edge_0 is the parent node of edge_1
    '''
    modify_tree = tmp_tree.copy("newick")
    if edge_0 != edge_1:
        new_node = Tree()
        new_node.add_child(name=new_add_node_name)
        detached_node = modify_tree&edge_1
        detached = detached_node.detach()
        inserted_node = modify_tree&edge_0
        inserted_node.add_child(new_node)
        new_node.add_child(detached_node)

    else:
        modify_tree.add_child(name=new_add_node_name)
        
    return modify_tree


def search_this_branch(tmp_tree, edge_0, edge_1, current_quartets, current_leave_node_name, queue, dic_for_leave_node_comb_name):
    modify_tree = get_modify_tree(tmp_tree, edge_0, edge_1, current_leave_node_name)
    modify_tree.resolve_polytomy(recursive=True)
    modify_tree.unroot()
    tmp_tree_score = judge_tree_score(modify_tree, current_quartets, current_leave_node_name, dic_for_leave_node_comb_name)

    dic = {}
    dic['tree'] = modify_tree
    dic['score'] = tmp_tree_score
    queue.put(dic)


def select_mask_node_pair(dl_predict, new_add_taxa):
    
    if new_add_taxa <= 9:
        return None

    mask_node_pair = []
    
    current_start = start_end_list[new_add_taxa][0]
    current_end = start_end_list[new_add_taxa][1]
    select_distribution = dl_predict[current_start:current_end+1]
    if np.max(select_distribution) < 0.90:
        return None
    else:
        x,y = nlargest_indices(select_distribution, int(max(10,0.01*len(select_distribution)))) 
        
    for i in range(0,len(x)):
        idx = x[i]
        topology_value = y[i]
        quartet_comb = comb_of_id[current_start+idx]
        
        if topology_value == 0:
            mask_node_pair.append((quartet_comb[0],quartet_comb[1]))
        if topology_value == 1:
            mask_node_pair.append((quartet_comb[0],quartet_comb[2]))
        if topology_value == 2:
            mask_node_pair.append((quartet_comb[1],quartet_comb[2]))
            
    return mask_node_pair


def mask_edge(tree,node1,node2,edge_list):
    # mask edge between node1 and node2
        
    if len(edge_list) <= 3:
        return edge_list

    ancestor_name = tree.get_common_ancestor(node1,node2).name
    remove_edge = []

    node = tree.search_nodes(name=node1)[0]
    while node:
        if node.name == ancestor_name:
            break

        edge_0 = node.up.name
        edge_1 = node.name

        if len(remove_edge) >= len(edge_list) - 3:
            break
        remove_edge.append((edge_0, edge_1))
        node = node.up

    node = tree.search_nodes(name=node2)[0]
    while node:
        if node.name == ancestor_name:
            break

        edge_0 = node.up.name
        edge_1 = node.name

        if len(remove_edge) >= len(edge_list) - 3:
            break
        remove_edge.append((edge_0, edge_1))
        node = node.up  

    for ele in remove_edge:
        if ele in edge_list:
            edge_list.remove(ele)
    
    return edge_list


def gen_phylogenetic_tree(current_quartets, beam_size):
    '''
        search the phylogenetic tree having highest score
        idx: the name of numpy file 
    '''
    current_leave_node_name = [chr(ord(u'\u4e00')+i) for i in range(0, taxa_num)]
    
    candidate_tree_beam = []

    quartet_id = leave_node_comb_name[0]

    for _label in [0, 1, 2]:
        if _label == 0:
            label_id = "".join([quartet_id[0], quartet_id[1], quartet_id[2], quartet_id[3]])
        elif _label == 1:
            label_id = "".join([quartet_id[0], quartet_id[2], quartet_id[1], quartet_id[3]])
        elif _label == 2:
            label_id = "".join([quartet_id[0], quartet_id[3], quartet_id[1], quartet_id[2]])
        
        _tree = tree_from_quartet(label_id)
        _tree.unroot()

        _tree_score = current_quartets[0, _label]

        tmp_tree_dict = {'Tree':_tree, 'tree_score':_tree_score}
        candidate_tree_beam.append(tmp_tree_dict)

        candidate_tree_beam.sort(key=lambda k: -k['tree_score'])

    idx_for_internal_node_name_pool = 0

    current_tree_score_beam = []
    optim_tree_beam = []

    #in the start point set beam size equal to 3
    for i in range(0, 3):
        current_tree_score_beam.append(candidate_tree_beam[i]['tree_score'])
        optim_tree_beam.append(candidate_tree_beam[i]['Tree'])
    
    for i in range(4, len(current_leave_node_name)):
        candidate_tree_beam = []

        for j in range(0, len(optim_tree_beam)):
            ele = optim_tree_beam[j]

            idx_for_this_iter = 0

            edge_0_list = []
            edge_1_list = []

            if ele == None:
                continue
            optim_tree = ele.copy("newick")

            for node in optim_tree.iter_descendants():
                tmp_tree = optim_tree.copy("newick")
                edge_0 = node.up.name
                edge_1 = node.name
                if edge_0 == '' or edge_1 == '':
                    continue

                else:
                    edge_0_list.append(edge_0)
                    edge_1_list.append(edge_1)

            queue = multiprocessing.Manager().Queue()
            para_list = [(tmp_tree, edge_0_list[k], edge_1_list[k], current_quartets, current_leave_node_name[i], queue, dic_for_leave_node_comb_name) for k in range(0, len(edge_0_list))]

            process_num = min(64, 2*i+6)
            pool = Pool(process_num)
            pool.starmap(search_this_branch, para_list)
            pool.close()
            pool.join()

            candidate_tree_score = []
            candidate_tree = []

            while not queue.empty():
                tmp_dic = queue.get()
                
                tmp_tree_dict = {'Tree':tmp_dic['tree'], 'tree_score':tmp_dic['score']+current_tree_score_beam[j]}
                candidate_tree_beam.append(tmp_tree_dict)

        candidate_tree_beam.sort(key=lambda k: -k['tree_score'])
        candidate_tree_beam = candidate_tree_beam[0:beam_size]
        
        optim_tree_beam = []
        current_tree_score_beam = []
        for ele in candidate_tree_beam:
            crt_tree = ele['Tree'].copy("newick") 
            for node in crt_tree.traverse("preorder"):
                if node.name == '':
                    node.name = str(internal_node_name_pool[idx_for_internal_node_name_pool])
                    idx_for_internal_node_name_pool += 1
            optim_tree_beam.append(crt_tree)
            crt_tree_score = ele['tree_score']
            current_tree_score_beam.append(crt_tree_score)

    return optim_tree_beam[0].write(format=9)


def gen_phylogenetic_tree_2(current_quartets, beam_size):
    '''
        search the phylogenetic tree having highest score
        idx: the name of numpy file 
    '''
    current_leave_node_name = [chr(ord(u'\u4e00')+i) for i in range(0, taxa_num)]
    
    candidate_tree_beam = []

    quartet_id = leave_node_comb_name[0]

    for _label in [0, 1, 2]:
        if _label == 0:
            label_id = "".join([quartet_id[0], quartet_id[1], quartet_id[2], quartet_id[3]])
        elif _label == 1:
            label_id = "".join([quartet_id[0], quartet_id[2], quartet_id[1], quartet_id[3]])
        elif _label == 2:
            label_id = "".join([quartet_id[0], quartet_id[3], quartet_id[1], quartet_id[2]])
        
        _tree = tree_from_quartet(label_id)
        _tree.unroot()

        _tree_score = current_quartets[0, _label]

        tmp_tree_dict = {'Tree':_tree, 'tree_score':_tree_score}
        candidate_tree_beam.append(tmp_tree_dict)

        candidate_tree_beam.sort(key=lambda k: -k['tree_score'])

    idx_for_internal_node_name_pool = 0

    current_tree_score_beam = []
    optim_tree_beam = []

    #in the start point set beam size equal to 3
    for i in range(0, 3):
        current_tree_score_beam.append(candidate_tree_beam[i]['tree_score'])
        optim_tree_beam.append(candidate_tree_beam[i]['Tree'])
    
    for i in range(4, len(current_leave_node_name)):
        candidate_tree_beam = []

        for j in range(0, len(optim_tree_beam)):
            ele = optim_tree_beam[j]

            idx_for_this_iter = 0

            edge_0_list = []
            edge_1_list = []

            if ele == None:
                continue
            optim_tree = ele.copy("newick")

            for node in optim_tree.iter_descendants():
                tmp_tree = optim_tree.copy("newick")
                edge_0 = node.up.name
                edge_1 = node.name
                if edge_0 == '' or edge_1 == '':
                    continue

                else:
                    edge_0_list.append(edge_0)
                    edge_1_list.append(edge_1)

            edge_list = [(edge_0_list[i], edge_1_list[i]) for i in range(0, len(edge_0_list))]
            mask_node_pairs = select_mask_node_pair(current_quartets, i)

            if mask_node_pairs != None:

                mask_node_pairs = list(set(mask_node_pairs))
                for node_pairs in mask_node_pairs:
                    node1 = chr(ord(u'\u4e00')+node_pairs[0])
                    node2 = chr(ord(u'\u4e00')+node_pairs[1])
                    
                    edge_list = mask_edge(ele.copy("deepcopy"),node1,node2,edge_list)
                    if len(edge_list) <= 3:
                        break

            edge_0_list = [ele[0] for ele in edge_list]
            edge_1_list = [ele[1] for ele in edge_list]

            queue = multiprocessing.Manager().Queue()
            para_list = [(tmp_tree, edge_0_list[k], edge_1_list[k], current_quartets, current_leave_node_name[i], queue, dic_for_leave_node_comb_name) for k in range(0, len(edge_0_list))]

            process_num = min(64, len(edge_0_list))
            pool = Pool(process_num)
            pool.starmap(search_this_branch, para_list)
            pool.close()
            pool.join()

            candidate_tree_score = []
            candidate_tree = []

            while not queue.empty():
                tmp_dic = queue.get()
                
                tmp_tree_dict = {'Tree':tmp_dic['tree'], 'tree_score':tmp_dic['score']+current_tree_score_beam[j]}
                candidate_tree_beam.append(tmp_tree_dict)

        candidate_tree_beam.sort(key=lambda k: -k['tree_score'])
        candidate_tree_beam = candidate_tree_beam[0:beam_size]
        
        optim_tree_beam = []
        current_tree_score_beam = []
        for ele in candidate_tree_beam:
            crt_tree = ele['Tree'].copy("newick") 
            for node in crt_tree.traverse("preorder"):
                if node.name == '':
                    node.name = str(internal_node_name_pool[idx_for_internal_node_name_pool])
                    idx_for_internal_node_name_pool += 1
            optim_tree_beam.append(crt_tree)
            crt_tree_score = ele['tree_score']
            current_tree_score_beam.append(crt_tree_score)

    return optim_tree_beam[0].write(format=9)


def transform_str(str_a, taxa_name):
    str_b = ''
    id_set = [chr(ord(u'\u4e00')+i) for i in range(0, taxa_num)]

    for i in range(0, len(str_a)):
        if str_a[i] in id_set:
            str_b += taxa_name[ord(str_a[i])-ord(u'\u4e00')]
        else:
            str_b += str_a[i]
    
    return str_b


def cmp(a, b):
    if int(b) > int(a):
        return -1
    if int(a) < int(b):
        return 1
    return 0


def get_numpy(aln_file):
    '''
        current version only supports the total length of msa less than 10K
    '''
    aln = open(aln_file)
    dic = {'A':'0','T':'1','C':'2','G':'3','-':'4', 'N':'4'}

    # for masking other unknown bases
    other_base = ['R', 'Y', 'K', 'M', 'U', 'S', 'W', 'B', 'D', 'H', 'V', 'X']
    for ele in other_base:
        dic[ele] = '4'

    matrix_out=[]
    fasta_dic={}
    for line in aln:
        if line[0]==">":
            header=line[1:].rstrip('\n').strip()
            fasta_dic[header]=[]
        elif line[0].isalpha() or line[0]=='-':
            for base, num in dic.items():
                line=line[:].rstrip('\n').strip().replace(base,num)
            line=list(line)
            line=[int(n) for n in line]
            fasta_dic[header] += line+[4]*(14000-len(line))

    taxa_block=[]
    for taxa in sorted(list(fasta_dic.keys()), key=functools.cmp_to_key(cmp)):
        taxa_block.append(fasta_dic[taxa.strip()])
    fasta_dic={}
    matrix_out.append(taxa_block)

    return np.array(matrix_out)


# def gen_quartet_seq(str, id, start_point):
#     index = np.array(str)
#     end_point = start_point + 1200
#     current_quartet = org_seq[0,index[:],start_point:end_point]
#     quartet_seq[id, comb_of_id.index(str)] = current_quartet


# def gen_quartet_seq_2(str, id, start_point):
#     index = np.array(str)
#     end_point = start_point + 240
#     current_quartet = org_seq[0,index[:],start_point:end_point]
#     quartet_seq[id, comb_of_id.index(str)] = current_quartet


def get_dl_model_1200():
    '''
    get the definition of dl model 1200
    this model aims to solve the default case 
    which are length larger than 1200 
    '''
    conv_x=[4,1,1,1,1,1,1,1]
    conv_y=[1,2,2,2,2,2,2,2]
    pool=[1,4,4,4,2,2,2,1]
    filter_s=[1024,1024,128,128,128,128,128,128]

    visible = layers.Input(shape=(4,1200,1))
    x = visible
        
    for l in list(range(0,8)):
        x = layers.ZeroPadding2D(padding=((0, 0), (0,conv_y[l]-1)))(x)        
        x = layers.Conv2D(filters=filter_s[l], kernel_size=(conv_x[l], conv_y[l]), strides=1, activation='relu')(x)
        x = layers.BatchNormalization()(x)
        x = layers.Dropout(rate=0.2)(x)
        x = layers.AveragePooling2D(pool_size=(1,pool[l]))(x)
        
    flat = layers.Flatten()(x)

    y = tf.keras.layers.Reshape((4,1200))(visible)
    y = tf.keras.layers.Bidirectional(tf.keras.layers.LSTM(128,return_sequences=True))(y)
    y = tf.keras.layers.Bidirectional(tf.keras.layers.LSTM(128,return_sequences=True))(y)
    y = tf.keras.layers.Bidirectional(tf.keras.layers.LSTM(128))(y)
    flat = tf.keras.layers.concatenate([flat, y],axis=-1)

    hidden1 = layers.Dense(1024,activation='relu')(flat)
    drop1 = layers.Dropout(rate=0.2)(hidden1)
    output = layers.Dense(3, activation='softmax')(drop1)
    model = tf.keras.Model(inputs=visible, outputs=output)

    return model


def get_dl_model_240():
    '''
    get the definition of dl model 240
    this model aims to solve the short length case 
    which are length larger than 240 
    '''
    conv_x=[4,1,1,1,1,1,1,1]
    conv_y=[1,2,2,2,2,2,2,2]
    pool=[1,2,2,2,2,2,2,2]
    filter_s=[1024,1024,128,128,128,128,128,128]

    visible = layers.Input(shape=(4,240,1))
    x = visible

    for l in list(range(0,8)):
        x = layers.ZeroPadding2D(padding=((0, 0), (0,conv_y[l]-1)))(x)        
        x = layers.Conv2D(filters=filter_s[l], kernel_size=(conv_x[l], conv_y[l]), strides=1, activation='relu')(x)
        x = layers.BatchNormalization()(x)
        x = layers.Dropout(rate=0.2)(x)
        x = layers.AveragePooling2D(pool_size=(1,pool[l]))(x)
        #print(x.shape)

    flat = layers.Flatten()(x)

    y = tf.keras.layers.Reshape((4,240))(visible)
    y = tf.keras.layers.Bidirectional(tf.keras.layers.LSTM(128,return_sequences=True))(y)
    y = tf.keras.layers.Bidirectional(tf.keras.layers.LSTM(128,return_sequences=True))(y)
    y = tf.keras.layers.Bidirectional(tf.keras.layers.LSTM(128))(y)
    flat = tf.keras.layers.concatenate([flat, y],axis=-1)

    hidden1 = layers.Dense(1024,activation='relu')(flat)
    drop1 = layers.Dropout(rate=0.2)(hidden1)
    output = layers.Dense(3, activation='softmax')(drop1)
    model = tf.keras.Model(inputs=visible, outputs=output)

    return model


def fill_dl_predict_each_slide_window(len_idx_1, len_idx_2):
    start_pos = 0
    iters = len(comb_of_id) // 50000
    for i in range(0, iters):
        batch_seq = np.zeros((50000, 4, 1200))

        for j in range(0, len(batch_seq)):
            idx = np.array(comb_of_id[i*50000+j])
            batch_seq[j] = org_seq[0, idx[:], len_idx_1:len_idx_2]

        test_seq = tf.expand_dims(batch_seq.astype(np.float32), axis=-1)
        predicted = dl_model.predict(x=test_seq)

        for j in range(0, len(batch_seq)):
            dl_predict[i*50000+j,:] += predicted[j,:]
        #dl_predict[i*50000:i*50000+len(batch_seq),:] += predicted[:,:]

        start_pos += 50000

    last_batch_size = len(comb_of_id) % 50000
    batch_seq = np.zeros((last_batch_size, 4, 1200))

    for j in range(0, len(batch_seq)):
        idx = np.array(comb_of_id[iters*50000+j])
        batch_seq[j] = org_seq[0, idx[:], len_idx_1:len_idx_2]

    test_seq = tf.expand_dims(batch_seq.astype(np.float32), axis=-1)
    predicted = dl_model.predict(x=test_seq)

    for j in range(0, len(batch_seq)):
        dl_predict[iters*50000+j,:] += predicted[j,:]
    #dl_predict[iters*50000:iters*50000+len(batch_seq),:] += predicted[:,:]


def fill_dl_predict_each_slide_window_2(len_idx_1, len_idx_2):
    start_pos = 0
    iters = len(comb_of_id) // 50000
    for i in range(0, iters):
        batch_seq = np.zeros((50000, 4, 240))

        for j in range(0, len(batch_seq)):
            idx = np.array(comb_of_id[i*50000+j])
            batch_seq[j] = org_seq[0, idx[:], len_idx_1:len_idx_2]

        test_seq = tf.expand_dims(batch_seq.astype(np.float32), axis=-1)
        predicted = dl_model.predict(x=test_seq)

        for j in range(0, len(batch_seq)):
            dl_predict[i*50000+j,:] += predicted[j,:]
        #dl_predict[i*50000:i*50000+len(batch_seq),:] += predicted[:,:]

        start_pos += 50000

    last_batch_size = len(comb_of_id) % 50000
    batch_seq = np.zeros((last_batch_size, 4, 240))

    for j in range(0, len(batch_seq)):
        idx = np.array(comb_of_id[iters*50000+j])
        batch_seq[j] = org_seq[0, idx[:], len_idx_1:len_idx_2]

    test_seq = tf.expand_dims(batch_seq.astype(np.float32), axis=-1)
    predicted = dl_model.predict(x=test_seq)

    for j in range(0, len(batch_seq)):
        dl_predict[iters*50000+j,:] += predicted[j,:]
    #dl_predict[iters*50000:iters*50000+len(batch_seq),:] += predicted[:,:]


def fill_dl_predict(window_number):
    step = (len_of_msa - 1200) // window_number
    start_idx = 0
    for i in range(0, window_number):
        end_idx = start_idx + 1200
        fill_dl_predict_each_slide_window(start_idx, end_idx)
        start_idx += step


def fill_dl_predict_2(window_number):
    if len_of_msa > 240:
        step = (len_of_msa - 240) // window_number
        start_idx = 0
        for i in range(0, window_number):
            end_idx = start_idx + 240
            fill_dl_predict_each_slide_window_2(start_idx, end_idx)
            start_idx += step

    else:
        step = (250 - 240) // window_number
        start_idx = 0
        for i in range(0, window_number):
            end_idx = start_idx + 240
            fill_dl_predict_each_slide_window_2(start_idx, end_idx)
            start_idx += step


if __name__ == '__main__':
    parser = argparse.ArgumentParser('get_msa_dir')
    p_input = parser.add_argument_group("INPUT")
    p_input.add_argument("--msa_dir", action="store", type=str, required=True)
    p_input.add_argument("--save_prefix", action="store", type=str, required=True)
    p_input.add_argument("--beam_size", action="store", type=str, default='3', required=False)
    p_input.add_argument("--sequence_type", action="store", type=str, default='standard', required=False)
    p_input.add_argument("--branch_model", action="store", type=str, default='gamma', required=False)
    p_input.add_argument("--window_coverage", action="store", type=str, default='1', required=False)


    args = parser.parse_args()
    msa_dir = args.msa_dir
    save_prefix = args.save_prefix
    beam_size = args.beam_size
    sequence_type = args.sequence_type
    branch_model = args.branch_model
    window_coverage = args.window_coverage

    flag = 0
    support_format = ['.fas', '.phy', '.fasta', 'phylip']
    bio_format = ['fasta', 'phylip', 'fasta', 'phylip']

    taxa_name = {}

    for i in range(0, len(support_format)):
        ele = support_format[i]
        if msa_dir.endswith(ele):
            flag = 1
            try:
                alignment = AlignIO.read(open(msa_dir), bio_format[i])

                len_of_msa = len(alignment[0].seq)
                taxa_num = len(alignment)

                save_alignment = save_prefix + '_fusang.fas'
                with open(save_alignment,'w') as f:
                    for record in alignment:
                        taxa_name[len(taxa_name)] = record.id
                        f.write('>'+str(len(taxa_name)-1)+'\n')
                        f.write(str(record.seq)+'\n')
            except: 
                print('Something wrong about your msa file, please check your msa file')
            break

    if flag == 0:
        print('we do not support this format of msa')
        sys.exit(1)

    #logging.warning(taxa_num)

    start_end_list = [None, None, None]

    end = -1
    for i in range(3, 100):
        start = end + 1
        end = start + int(comb_math(i,3)) - 1
        start_end_list.append((start,end))

    id_for_taxa = [i for i in range(0, taxa_num)]
    comb_of_id = list(combinations(id_for_taxa, 4))
    comb_of_id.sort(key=lambda ele: ele[-1])

    leave_node_id = [i for i in range(0, taxa_num)]
    leave_node_name = [chr(ord(u'\u4e00')+i) for i in range(0, taxa_num)]
    leave_node_comb_id = comb_of_id

    leave_node_comb_name = []

    dic_for_leave_node_comb_name = {}

    for ele in leave_node_comb_id:
        term = [chr(ord(u'\u4e00')+id) for id in ele]
        dic_for_leave_node_comb_name["".join(term)] = len(dic_for_leave_node_comb_name)
        leave_node_comb_name.append("".join(term))
    
    internal_node_name_pool = ['internal_node_' + str(i) for i in range(3, 3000)]

    fusang_msa_dir = save_prefix + '_fusang.fas'
    org_seq = get_numpy(fusang_msa_dir)
    os.remove(fusang_msa_dir)

    window_number = 1

    if len_of_msa <= 1210:
        dl_model = get_dl_model_240()
        if sequence_type == 'standard' and branch_model == 'gamma':
            dl_model.load_weights(filepath='./dl_model/len_240/S1G/best_weights_clas').expect_partial()
        if sequence_type == 'standard' and branch_model == 'uniform':
            dl_model.load_weights(filepath='./dl_model/len_240/S1U/best_weights_clas').expect_partial()
        if sequence_type == 'coding' and branch_model == 'gamma':
            dl_model.load_weights(filepath='./dl_model/len_240/C1G/best_weights_clas').expect_partial()
        if sequence_type == 'coding' and branch_model == 'uniform':
            dl_model.load_weights(filepath='./dl_model/len_240/C1U/best_weights_clas').expect_partial()
        if sequence_type == 'noncoding' and branch_model == 'gamma':
            dl_model.load_weights(filepath='./dl_model/len_240/N1G/best_weights_clas').expect_partial()
        if sequence_type == 'noncoding' and branch_model == 'uniform':
            dl_model.load_weights(filepath='./dl_model/len_240/N1U/best_weights_clas').expect_partial()

        window_number = int(len_of_msa * float(window_coverage) // 240 + 1) 

    elif len_of_msa > 1210:
        dl_model = get_dl_model_1200()
        if sequence_type == 'standard' and branch_model == 'gamma':
            dl_model.load_weights(filepath='./dl_model/len_1200/S2G/best_weights_clas').expect_partial()
        if sequence_type == 'standard' and branch_model == 'uniform':
            dl_model.load_weights(filepath='./dl_model/len_1200/S2U/best_weights_clas').expect_partial()
        if sequence_type == 'coding' and branch_model == 'gamma':
            dl_model.load_weights(filepath='./dl_model/len_1200/C2G/best_weights_clas').expect_partial()
        if sequence_type == 'coding' and branch_model == 'uniform':
            dl_model.load_weights(filepath='./dl_model/len_1200/C2U/best_weights_clas').expect_partial()
        if sequence_type == 'noncoding' and branch_model == 'gamma':
            dl_model.load_weights(filepath='./dl_model/len_1200/N2G/best_weights_clas').expect_partial()
        if sequence_type == 'noncoding' and branch_model == 'uniform':
            dl_model.load_weights(filepath='./dl_model/len_1200/N2U/best_weights_clas').expect_partial()

        window_number = int(len_of_msa * float(window_coverage) // 1200 + 1) 


    dl_predict = np.zeros((len(comb_of_id), 3))

    if len_of_msa > 1210:
        fill_dl_predict(window_number)
    elif len_of_msa <= 1210:
        fill_dl_predict_2(window_number)
    else:
        print('current version of fusang do not support this length of MSA')
        sys.exit(1)

    dl_predict /= window_number

    #np.save('./dl_predict.npy', dl_predict)
    #dl_predict = np.load('./dl_predict.npy')

    if not os.path.exists('./dl_output/'):
        os.mkdir('./dl_output/')

    if taxa_num > 10:
        searched_tree = transform_str(gen_phylogenetic_tree_2(dl_predict, int(beam_size)), taxa_name)
    else:
        searched_tree = transform_str(gen_phylogenetic_tree(dl_predict, int(beam_size)), taxa_name)

    build_log = open('./dl_output/{}.txt'.format(save_prefix), 'a')
    build_log.write(searched_tree)
    build_log.close()

