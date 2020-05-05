import numpy as np
import pickle
import spec_tool
from spec_tool import SAMPLE_CHEM
from sklearn.cluster import DBSCAN
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import multiprocessing
from multiprocessing.managers import BaseManager

total_cpu = multiprocessing.cpu_count()  #在此设置本程序中可以调用的cpu总数
path = 'data'  #数据所在文件夹,注意末尾不用带 /


class Multi(object):
    threshold = 0.5  #阈值设定为0.5

    def __init__(self):
        self.start = 0  # 用来存储开始的下标
        self.process = ''  # 用来存储进程
        self.result = []  # 用来存储运算结果


def dealpeaklist(sample):
    if len(sample.Peak) == 0:
        sample.Peak = []
        return
    listone = sample.Peak.split(" ")
    peaklsit = []
    for i in listone:
        j = i.split(":")
        peaklsit.append((float(j[0]), float(j[1])))
    sample.Peak = peaklsit


def load_original_file():
    '''加载原始文件,返回由多个参数组成的字典'''
    global path
    index = 0
    blank, area_matrix = [], []
    with open(path + "/load_ontology.txt", 'rt') as txtfile:
        ontology_list = txtfile.read().split("\n")
    chemnum_list = [i for i in range(len(ontology_list))]
    with open(path + "/load_area.txt", "rt") as txtfile:
        for lines in txtfile:
            lines = lines.replace("\n", '')
            if index < 4:
                if index == 0:
                    speciesnum_list = lines.split("\t")[1:]
                elif index == 1:
                    species_list = lines.split("\t")[1:]
                elif index == 2:
                    scientific_list = lines.split("\t")[1:]
                elif index == 3:
                    clade_list = lines.split("\t")[1:]
                index += 1
            else:
                lines = lines.split("\t")
                blank.append(lines[0])
                area_matrix.append(lines[1:])
    blank = np.array(blank, dtype='float32')
    area_matrix = np.array(area_matrix, dtype='int32')

    sample_list = []
    with open(path + "/load_peak.txt", 'rt') as txtfile:
        for lines in txtfile:
            lines = lines.replace("\n", '').split("\t")
            sample = SAMPLE_CHEM()
            i = 0
            for key in sample.__dict__:
                sample.__dict__[key] = lines[i]
                i += 1
            dealpeaklist(sample)
            sample_list.append(sample)

    for s in sample_list:
        if len(s.Peak) == 0:
            continue
        spec_tool.addinformation(s)

    total_dict = {
        "speciesnum_list": speciesnum_list,  # 物种序号
        "species_list": species_list,  # 中文物种名
        "scientific_list": scientific_list,  # 拉丁物种名
        "clade_list": clade_list,  # 物种的clade
        "area_matrix": area_matrix,  # 代谢物含量的矩阵
        "sample_list": sample_list,  # 化合物的峰
        "blank": blank,  #空白对照
        "ontology_list": ontology_list,  # 化合物的ontology
        "chemnum_list": chemnum_list  # 化合物序号
    }
    return total_dict


def delete_impurity(total_dict):
    '''删去杂质'''
    area_matrix = total_dict['area_matrix']
    sample_list = total_dict['sample_list']
    blank = total_dict['blank']
    max_list = np.max(area_matrix, axis=1)  #对每一行求最大值
    index = 0
    while index < len(area_matrix):
        if max_list[index] < blank[index] * 10 or 'Stearamide' in sample_list[
                index].NAME or 'Erucamide' in sample_list[
                    index].NAME or 'Dimethoate' in sample_list[index].NAME:
            area_matrix = np.delete(area_matrix, index, axis=0)
            blank = np.delete(blank, index, axis=0)
            max_list = np.delete(max_list, index, axis=0)
            del sample_list[index]
            del total_dict['ontology_list'][index]
            del total_dict['chemnum_list'][index]
            # 如果均值没有比空白高10倍,则认为是杂质,删除之
            # 要指定删除的名字: Stearamide,Erucamide
        else:
            index += 1
    del total_dict['speciesnum_list']
    del total_dict['species_list']
    del total_dict['blank']  # 删除不会用到的部分,减少内存压力


def my_DBSCAN(log_matrix):
    print("To my_DBSCAN")
    high_lines_clusters, low_lines_clusters = [], []
    for line in log_matrix:
        this_line_high, this_line_low = [], []
        X = np.array(line).reshape(-1, 1)
        label_list = list(DBSCAN(eps=0.5, min_samples=1).fit(X).labels_)
        num_of_label = len(set(label_list))  #总共分为了这么多类
        represent = []
        now_label = 0
        while now_label < num_of_label:
            for i in range(len(label_list)):
                if label_list[i] == now_label:
                    represent.append(line[i])
                    now_label += 1
                    break
        label_order = np.argsort(represent)[::-1]
        this_line_pointer = this_line_high
        for label in label_order:
            # 不取组内数量大于总体组数1/6的组,将含量特别高和含量特别低的那些分开两组
            if label_list.count(label) > len(log_matrix[0]) / 6:
                this_line_pointer = this_line_low
                continue
            cluster_list = []
            for index in range(len(label_list)):
                if label_list[index] == label:
                    cluster_list.append(index)  #此时index与对应物种序列号相等
            this_line_pointer.append(cluster_list)
        high_lines_clusters.append(this_line_high)
        low_lines_clusters.append(this_line_low)
    for index in range(len(low_lines_clusters)):
        low_lines_clusters[index] = low_lines_clusters[index][::-1]
        # 跟high_lines_clusters统一,都是将最极端的放在最前面
    return high_lines_clusters, low_lines_clusters


def chemistry_similarity(result, sample_list, start, limit):
    '''计算化合物之间的两两相似度,该函数被每个单一进程独立调用'''
    threshold = Multi.threshold
    similarity_matrix = []
    for i in range(start, start + limit):
        if i >= len(sample_list):
            break
        if len(sample_list[i].Peak) == 0:
            similarity_matrix.append([0 for _ in range(len(sample_list))])
            continue
        zero_list = [0 for _ in range(i + 1)]
        this_line = []
        for j in range(i + 1, len(sample_list)):
            if len(sample_list[j].Peak) == 0:
                this_line.append(0)
                continue
            score = spec_tool.mix_similarity(sample_list[i], sample_list[j])
            score = score if score >= threshold else 0
            this_line.append(score)
        zero_list.extend(this_line)
        similarity_matrix.append(zero_list)

    similarity_matrix = np.array(similarity_matrix, dtype='float32')
    for i in range(len(similarity_matrix)):
        for j in range(len(similarity_matrix[i])):
            if similarity_matrix[i][j] > 0.001:
                similarity_matrix[i][j] = (similarity_matrix[i][j] -
                                           threshold) * (1 / (1 - threshold))
    # 转为0-1
    result.append(similarity_matrix)


def parallel_computing(sample_list, distance_matrix, high_lines_clusters,
                       start_index, limit):
    '''使用并行运算来计算化合物之间的相似性,在该函数中分配资源,在chemistry_similarity中实现具体计算步骤'''
    process_list = []
    lock = multiprocessing.Lock()
    with multiprocessing.Manager() as manager:  # 使用该方法保证列表能在进程间传递
        global total_cpu
        for cpu in range(total_cpu - 1):
            if start_index >= len(sample_list):
                break
            m = Multi()
            m.start = start_index
            m.result = manager.list()
            m.process = multiprocessing.Process(target=chemistry_similarity,
                                                args=(m.result, sample_list,
                                                      m.start, limit))
            start_index += limit
            m.process.start()
            process_list.append(m)
        for m in process_list:
            m.process.join()
            m.result = m.result[0]
            lock.acquire()
            similarity_distance(m, distance_matrix,
                                high_lines_clusters[m.start:limit])
            del m.result  # 减轻内存压力
            lock.release()
    return distance_matrix


def DBSCAN_distance(high_lines_clusters, low_lines_clusters, species_num):
    print("To calculate_distance")
    zero_list = [0 for _ in range(species_num)]
    distance_matrix = [zero_list[::] for _ in range(species_num)]
    distance_matrix = np.array(distance_matrix, dtype='float32')
    for all_lines_clusters in high_lines_clusters, low_lines_clusters:
        for each_line in all_lines_clusters:
            # 先计算同一类的
            if len(each_line) == 0:
                continue
            affinity = 1 + 1 / len(each_line)
            for cell in each_line:
                affinity -= 1 / len(each_line)
                for i in range(len(cell) - 1):
                    for j in range(i + 1, len(cell)):
                        distance_matrix[cell[i]][cell[j]] += affinity
                        distance_matrix[cell[j]][cell[i]] += affinity
            # 再计算类与类之间的
            for i in range(len(each_line) - 1):
                for j in range(i + 1, len(each_line)):
                    affinity = 1 - 1 / len(each_line) * (j - 1)
                    for first in each_line[i]:
                        for second in each_line[j]:
                            distance_matrix[first][second] += affinity
                            distance_matrix[second][first] += affinity
    return distance_matrix


def similarity_distance(multiobj, distance_matrix, high_lines_clusters):
    '''把物质之间的相似度计算进入样本距离,仅对含量高的物质计算距离'''
    for si in range(len(multiobj.result) - 1):
        for sj in range(si + 1, len(multiobj.result)):
            score = multiobj.result[si][sj]
            if score < 0.01 or len(high_lines_clusters[si]) == 0:
                continue
            affinity_i = 1 + 1 / len(high_lines_clusters[si])
            for cluster_one in high_lines_clusters[si]:
                affinity_i -= 1 / len(high_lines_clusters[si])
                if len(high_lines_clusters[sj]) == 0:
                    continue
                affintiy_j = 1 + 1 / len(high_lines_clusters[sj])
                for cluster_two in high_lines_clusters[sj]:
                    affintiy_j -= 1 / len(high_lines_clusters[sj])
                    for index_one in cluster_one:
                        for index_two in cluster_two:
                            distance_matrix[index_one][
                                index_two] += score * affinity_i * affintiy_j
                            distance_matrix[index_two][
                                index_one] += score * affinity_i * affintiy_j


def output_distance_matrix(distance_matrix):
    global path
    with open(path + "/distance_matrix.txt", 'wt') as txtfile:
        for line in distance_matrix:
            string = ''
            for num in line:
                string += str(num) + '\t'
            string += '\n'
            txtfile.write(string)


def find_max_affinity(distance_matrix):
    x, y, max_affinity = 0, 0, 0
    for i in range(len(distance_matrix) - 1):
        for j in range(i + 1, len(distance_matrix)):
            if distance_matrix[i][j] > max_affinity:
                x, y = i, j
                max_affinity = distance_matrix[i][j]
    return x, y, max_affinity  # 注意这里x比y小


def find_min_distance(new_distance_matrix):
    x, y, min_distance = 0, 0, float("inf")
    for i in range(len(new_distance_matrix) - 1):
        for j in range(i + 1, len(new_distance_matrix)):
            if new_distance_matrix[i][j] < min_distance:
                x, y = i, j
                min_distance = new_distance_matrix[i][j]
    return x, y, min_distance  # 注意这里x比y小


def yield_new_line(distance_matrix, x, y):
    new_line = []
    for i in range(len(distance_matrix)):
        if distance_matrix[x][i] < distance_matrix[y][i]:
            new_line.append(distance_matrix[y][i])
        else:
            new_line.append(distance_matrix[x][i])
    return new_line


def hierarchial_clustering(distance_matrix, species_list, area_matrix):
    print("To heierarchial_clustering")
    clustering_matrix = []
    species_length, original_length = len(species_list), len(species_list)
    number_list = [i for i in range(species_length)]
    _, __, for_distance = find_max_affinity(distance_matrix)
    area_matrix = list(np.array(area_matrix, dtype='float32').T)
    # 这里要转置一下,因为原来横行是每个代谢物在不同样本中的含量,现在要改为每个样本中不同代谢物的含量
    while len(distance_matrix) > 1 and distance_matrix.max(
    ) > distance_matrix.min():
        x, y, max_affinity = find_max_affinity(distance_matrix)
        x_name = '' if number_list[x] >= original_length else species_list[
            number_list[x]]
        y_name = '' if number_list[y] >= original_length else species_list[
            number_list[y]]
        xy_distance = for_distance - max_affinity + 1
        # 进化枝之间的距离定义为跟最高亲和度的差值+1
        new_line = yield_new_line(distance_matrix, x, y)
        distance_matrix = np.r_[distance_matrix, [new_line]]
        new_line.extend([0.0])
        distance_matrix = np.c_[distance_matrix, new_line]
        distance_matrix = np.delete(distance_matrix, y, axis=0)
        distance_matrix = np.delete(distance_matrix, x, axis=0)
        distance_matrix = np.delete(distance_matrix, y, axis=1)
        distance_matrix = np.delete(distance_matrix, x, axis=1)
        clustering_matrix.append(
            [number_list[x], number_list[y], xy_distance, x_name, y_name])
        del number_list[y], number_list[x]
        number_list.append(species_length)
        species_length += 1

        area_new_line = []
        for k in x, y:
            if isinstance(area_matrix[k][0], list) or isinstance(
                    area_matrix[k][0], np.ndarray):
                area_new_line.extend(area_matrix[k])
            else:
                area_new_line.append(area_matrix[k])
        area_matrix.append(area_new_line)

    if len(distance_matrix) > 1:
        new_area_matrix = []
        for index in range(len(distance_matrix)):
            current_line = area_matrix[number_list[index]]
            if isinstance(current_line[0], list) or isinstance(
                    current_line[0], np.ndarray):
                for index in range(1, len(current_line)):
                    current_line[0] += current_line[1]
                new_area_matrix.append(current_line[0] / len(current_line))
            else:
                new_area_matrix.append(current_line)

        new_area_matrix = np.array(new_area_matrix, dtype='float32')
        scaler = StandardScaler()  # 归一化处理
        new_area_matrix = scaler.fit(new_area_matrix).transform(
            new_area_matrix)
        pca = PCA(n_components=len(new_area_matrix) - 1)  # 主成分分析
        pca.fit(new_area_matrix)
        new_area_matrix = pca.transform(new_area_matrix)
        while len(new_area_matrix) > 1:
            inf_list = [float('inf') for _ in range(len(new_area_matrix))]
            new_distance_matrix = [
                inf_list[::] for _ in range(len(new_area_matrix))
            ]
            for i in range(len(new_area_matrix)):
                for j in range(len(new_area_matrix)):
                    new_distance_matrix[i][j] = np.linalg.norm(
                        new_area_matrix[i] - new_area_matrix[j])
            x, y, xy_distance = find_min_distance(new_distance_matrix)
            x_name = '' if number_list[x] > original_length else species_list[
                number_list[x]]
            y_name = '' if number_list[y] > original_length else species_list[
                number_list[y]]
            area_new_line = (new_area_matrix[x] + new_area_matrix[y]) / 2
            new_area_matrix = np.r_[new_area_matrix, [area_new_line]]
            new_area_matrix = np.delete(new_area_matrix, y, axis=0)
            new_area_matrix = np.delete(new_area_matrix, x, axis=0)
            clustering_matrix.append(
                [number_list[x], number_list[y], xy_distance, x_name, y_name])
            del number_list[y], number_list[x]
            number_list.append(species_length)
            species_length += 1

    return clustering_matrix


def choose_data(total_dict):
    '''可以指定某几个ontology进行代谢物的进化关系分析'''
    std = ['fla', 'antho']  # 要指定的ontology名的词根
    olist = total_dict['ontology_list']
    indexlist = [i for i in range(len(olist)) for s in std if s in olist[i]]
    new_samlist, new_arealist = [], []
    for i in indexlist:
        new_arealist.append(total_dict['area_matrix'][i])
        new_samlist.append(total_dict['sample_list'][i])
    total_dict['area_matrix'] = np.array(new_arealist)
    total_dict['sample_list'] = new_samlist
    return total_dict


def main():
    global path
    # total_dict = load_original_file()
    # with open(path+"/sample_dict.pkl", 'wb') as pklfile:
    #     pickle.dump(total_dict, pklfile)
    # return
    with open(path + "/sample_dict.pkl", "rb") as pklfile:
        total_dict = pickle.load(pklfile)
    delete_impurity(total_dict)

    # 取前300个样品,前80个物种来做测试
    # total_dict['sample_list'] = total_dict['sample_list'][:300]
    # total_dict['area_matrix'] = total_dict['area_matrix'][:300, :80]
    # total_dict['species_list'] = total_dict['species_list'][:80]

    # choose_data(total_dict)
    high_lines_clusters, low_lines_clusters = my_DBSCAN(
        np.log10(total_dict['area_matrix'] + 1))
    distance_matrix = DBSCAN_distance(high_lines_clusters, low_lines_clusters,
                                      len(total_dict['scientific_list']))

    global total_cpu
    start_index = 0
    limit = 100
    while start_index < len(total_dict['sample_list']):
        distance_matrix = parallel_computing(total_dict['sample_list'],
                                             distance_matrix,
                                             high_lines_clusters, start_index,
                                             limit)
        start_index += limit * (total_cpu - 1)
        print(start_index)
        ###################这里保存半成品,便于下次继续进行
        with open(path + "/half_result.pkl", 'wb') as pklfile:
            half_dict = {
                'distance_matrix': distance_matrix,
                'start_index': start_index,
                'high_lines_clusters': high_lines_clusters
            }
            pickle.dump(half_dict, pklfile)


def continue_main():
    global path
    global total_cpu
    with open("work/half_result.pkl", 'rb') as pklfile:
        half_dict = pickle.load(pklfile)  #接着上次的结果来做
    start_index = half_dict['start_index']
    distance_matrix = half_dict['distance_matrix']
    high_lines_clusters = half_dict['high_lines_clusters']
    with open(path + '/sample_dict.pkl', 'rb') as pklfile:
        total_dict = pickle.load(pklfile)

    limit = 300
    while start_index < len(total_dict['sample_list']):
        distance_matrix = parallel_computing(total_dict['sample_list'],
                                             distance_matrix,
                                             high_lines_clusters, start_index,
                                             limit)
        start_index += limit * (total_cpu - 1)
        print(start_index)
        ###################这里保存半成品,便于下次继续进行
        with open("work/half_result.pkl", 'wb') as pklfile:
            half_dict = {
                'distance_matrix': distance_matrix,
                'start_index': start_index,
                'high_lines_clusters': high_lines_clusters
            }
            pickle.dump(half_dict, pklfile)

    output_distance_matrix(distance_matrix)
    clustering_matrix = hierarchial_clustering(distance_matrix,
                                               total_dict['scientific_list'],
                                               total_dict['area_matrix'])
    fasttree_string = spec_tool.transform_fasttree(
        clustering_matrix, len(total_dict['scientific_list']))
    with open(path + "/clustering_result.nwk", 'wt') as nwkfile:
        nwkfile.write(fasttree_string)

    # 以下输出的是clade,clade的建树是用来给后面打分用的
    clade_matrix = hierarchial_clustering(distance_matrix,
                                          total_dict['clade_list'],
                                          total_dict['area_matrix'])
    clade_string = spec_tool.transform_fasttree(clade_matrix,
                                                len(total_dict['clade_list']))
    with open(path + "/clade_result.nwk", 'wt') as nwkfile:
        nwkfile.write(clade_string)


if __name__ == '__main__':
    main()
    # continue_main()