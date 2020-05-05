from spec_tool import SAMPLE_CHEM
import spec_tool
import numpy as np
from spec_tool import between
import pickle
import json

max_category=6 # 设置递归最深层数,如果递归层数过深,网络会十分冗杂,不利于进行分子网络可视化

def find_by_known(knownobj, unknownlist, category):
    '''knownobj为一个已知物质,unknolist为未知物质的列表,根据已知物质来搜索相似性高的未知物质,以余弦相似度为标准'''
    threshold = 0.6
    result = []
    for un in unknownlist:
        if un.marker or len(un.Peak)==0:
            continue  # 表明该物质已知或没有离子碎片
        if compare_weight(un.Weight, knownobj.Weight): #差距超过阈值
            continue
        score = spec_tool.mix_similarity(knownobj, un)
        un.score = score
        un.category = category
        if score > threshold:
            result.append(un)
            un.marker = True
    return result


def findnext(node, unknownlist, category):
    global max_category
    if category>max_category:
        return []
    node.next = find_by_known(node, unknownlist, category)
    if len(node.next) > 0:
        for nextnode in node.next:
            findnext(nextnode, unknownlist, category + 1)  # 递归往下找


def depth_first(treenode, links):
    ''''对树进行遍历,将所有结果都添加到列表中'''
    if hasattr(treenode,"next") and len(treenode.next) > 0:
        for nextnode in treenode.next:
            d = {
                "source": treenode.Alignment_ID,
                "target": nextnode.Alignment_ID,
                "lineStyle": {
                    "width": nextnode.score
                }
            }
            links.append(d)
            depth_first(nextnode, links)


def main():
    path = 'data/'
    with open(path + 'sample_dict.pkl', 'rb') as pklfile:
        sample_dict = pickle.load(pklfile)
    sample_list = sample_dict['sample_list']
    del sample_dict

    for s in sample_list:
        s.marker = False
        spec_tool.addinformation(s)
    for i in range(len(sample_list)):
        if int(sample_list[i].Alignment_ID) == 60389:
            sample_list[i].marker = True
            sample_list[i].category = 0
            break
    findnext(sample_list[i], sample_list, 1)

    nodes, links, categories = [], [], []
    depth_first(sample_list[i], links)
    maxcategory = 0
    for i in range(len(sample_list)):
        if sample_list[i].marker:
            nodes.append({
                "name": sample_list[i].Alignment_ID,
                "category": sample_list[i].category
            })
            maxcategory = max(maxcategory, int(sample_list[i].category))
    for i in range(maxcategory + 1):
        categories.append({"name": str(i), "keyword": {}, "base": str(i)})
    jsondict = {
        "type": "force",
        "categories": categories,
        "nodes": nodes,
        "links": links
    }
    with open(path + "Nobiletin_network.json", 'wt') as jsonfile:
        jsonfile.write(json.dumps(jsondict))


def compare_weight(weight1, weight2):
    distance = abs(weight1 - weight2)
    smaller = min(weight1, weight2)
    if distance / smaller > 3:
        return True
    else:
        return False


if __name__ == "__main__":
    main()
