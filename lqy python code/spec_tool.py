import numpy as np

'''该python为工具包,其他程序中常用函数统一汇集在此处,便于调用;本python中共有三个类对象,三个类对象反应的是不同的数据来源,在实际使用过程中可当做同一对象进行计算'''

class STDARD_CHEM(object):
    def __init__(self):
        self.NAME = ''
        self.PRECURSORMZ = 0.0
        self.PRECURSORTYPE = ''
        self.FORMULA = ''
        self.Ontology = ''
        self.INCHIKEY = ''
        self.SMILES = ''
        self.RETENTIONTIME = ''
        self.IONMODE = ''
        self.INSTRUMENTTYPE = ''
        self.INSTRUMENT = ''
        self.COLLISIONENERGY = ''
        self.Comment = ''
        self.NumPeaks = None
        self.Peak = []  # 元组格式:(分子量,丰度)

    def __str__(self):
        string = ''
        for key in self.__dict__:
            string += key + ":" + str(self.__dict__[key]) + '\n'
        return string


class TEST_CHEM(object):
    def __init__(self):
        self.Alignment_ID = ''  # 注意有五位数,不够五位的左边填充了0
        self.Time = 0.0  # 吸附时间
        self.NAME = ''  # 特别注意,如果是Unknown的实际是''
        self.Weight = 0.0
        self.Peak = []  # 元组格式:(分子量,丰度)

    def __str__(self):
        string = ''
        for key in self.__dict__:
            string += key + ":" + str(self.__dict__[key]) + '\n'
        return string


class SAMPLE_CHEM(object):
    def __init__(self):
        self.Alignment_ID = ''
        self.Time = ''
        self.Weight = ''
        self.NAME = ""
        self.Reference = ''
        self.Peak = []

    def __str__(self):
        string = ''
        for key in self.__dict__:
            string += key + ":" + str(self.__dict__[key]) + '\n'
        return string


def count_basepeak(peaklist):
    '''用于计算碎片的基峰,传入参数为峰的列表,返回丰度最高的峰'''
    if (peaklist == [] or peaklist == None):
        return None
    basepeak = peaklist[0][1]
    for each_peak in peaklist:
        if (each_peak[1] > basepeak):
            basepeak = each_peak[1]
    return basepeak


def count_maxweight(peaklist):
    '''用于计算碎片中最大分子量,传入参数为峰的列表,返回最大分子量'''
    if (peaklist == [] or peaklist == None):
        return None
    maxweight = peaklist[0][0]
    for each_peak in peaklist:
        if (each_peak[0] > maxweight):
            maxweight = each_peak[0]
        return maxweight


def drawpeak(peaklist_origin):
    '''该方法可绘制质谱图,传入参数为峰的列表,返回值为空'''
    def goto(turtle, x, y):
        '''由于每次转换位置都要抬起和放下画笔很麻烦，特此重写方法，传入的turtle为对象'''
        turtle.up()
        turtle.goto(x, y)
        turtle.down()

    if len(peaklist_origin) == 0:  #列表不空才能画
        print("这个物质的峰碎片为空列表")
        return
    peaklist = peaklist_origin[:]  #复制一个新表,避免画图对原表改动
    if isinstance(peaklist[0][0], str):
        peaklist = str_to_num(peaklist)  #如果是字符型的不能进行计算,会报出类型错误,因此要先进行类型转换
    basepeak = count_basepeak(peaklist)

    #以下变量为绘图基本元素,统一列于此处便于修改----------------------------
    window_width = 1000
    window_height = 600  #窗口的宽度和高度
    start_X = -400
    start_Y = -200  #坐标轴零点位置
    axis_size = 3  #画坐标轴时画笔粗细，范围为1到10
    length_X = 800
    length_Y = 450  #XY轴的长度
    axis_distance = 20  #数字与坐标轴的宽度间隙
    Y_gap = int(length_Y / 6)  #Y轴上标的单位长度
    #-------------------------------------------------------------------

    import turtle as t
    t.setup(window_width, window_height)
    t.hideturtle()
    t.speed('fastest')  #取消动画效果
    goto(t, start_X, start_Y)  #前往原点

    t.pensize(axis_size)
    t.forward(length_X)  #绘制X轴
    goto(t, start_X + length_X, start_Y + axis_distance)
    t.write("Weight", font=('Times New Roman', '12', 'normal'))
    goto(t, start_X, start_Y)
    t.left(90)
    t.forward(length_Y)  #绘制Y轴

    for number in range(1, 7):
        goto(t, start_X - axis_distance, start_Y + number * Y_gap)
        t.write(str(number * 20), font=('Times New Roman', '10', 'normal'))

    peaklist = sort_by_weight(peaklist)
    max_X = peaklist[0][0] * 1.2
    t.pensize(1)
    number_down = 1  #表示X轴上写的两个分子量要不要错开,1表示不用,2表示错开
    for each_peak in peaklist:  #画每一个峰
        X_position = each_peak[0] / max_X * length_X + start_X
        goto(t, X_position, start_Y)
        peak_height = each_peak[1] / basepeak * (Y_gap * 5)
        t.forward(peak_height)  #画出了峰的那条线

        goto(t, X_position, start_Y + peak_height)
        relative_peak = round(each_peak[1] / basepeak * 100, 3)
        t.write(str(relative_peak))  #写下丰度的数字

        index = peaklist.index(each_peak)  #下面这里是要判断两个数字之间要不要错开
        if (index != 0):
            weight_gap = max_X / 16
            if (abs(peaklist[index - 1][0] - peaklist[index][0]) < weight_gap):
                number_down = 2 if number_down == 1 else 1
            else:
                number_down = 1
        goto(t, X_position, start_Y - axis_distance * number_down)
        t.write(str(round(each_peak[0], 3)))  #写下分子量的数字

    t.done()


def sort_by_weight(peaklist):
    '''根据分子量从小到大为峰的列表排序'''
    for i in range(len(peaklist) - 1):
        flag = True
        for j in range(i, len(peaklist) - 1):
            try:
                if (peaklist[j][0] > peaklist[j + 1][0]):
                    peaklist[j], peaklist[j + 1] = peaklist[j + 1], peaklist[j]
                    flag = False
            except TypeError:
                peaklist[j + 1] = (float(peaklist[j + 1][0]),
                                   float(peaklist[j + 1][1]))
                peaklist[j] = (float(peaklist[j][0]), float(peaklist[j][1]))
                if (peaklist[j][0] > peaklist[j + 1][0]):
                    peaklist[j], peaklist[j + 1] = peaklist[j + 1], peaklist[j]
                    flag = False
        if flag:
            break  #flag仍然为真说明没有经过交换,也就是已经有序了
    return peaklist




def str_to_num(peaklist):
    '''该方法可以将以字符串形式存储的峰转换为数值型,返回值为转换后的峰'''
    num_peaklist = []
    for i in range(len(peaklist)):
        num_peaklist.append((float(peaklist[i][0]), float(peaklist[i][1])))
    return num_peaklist


def between(subject, measure, size=0.1):
    return measure - size <= subject <= measure + size


def relative_intensity(peaklist):
    '''将绝对丰度转化为相对丰度'''
    basepeak = count_basepeak(peaklist)
    return [(p[0], p[1] / basepeak) for p in peaklist]


def cos_distance(testobj1, testobj2):
    '''计算两个峰列表之间的余弦距离,传入参数为两个对象,注意需要包含相对丰度列表,W列表和距离'''
    W1, W2 = insert_zero(testobj1, testobj2)
    numerator = W1.dot(W2.T)
    denom = testobj1.distance * testobj2.distance
    if denom < 1:
        return 0
    else:
        return numerator / denom


def W(peaklist):
    '''传入一个相对丰度的峰列表,返回计算好的W值列表'''
    result_list = []
    for p in peaklist:
        result = p[0]**2 * np.sqrt(p[1])
        result_list.append(result)
    return result_list


def combine_peak(peaklist, size=0.01):
    '''将分子量在size之内的峰合并'''
    i = 0
    newlist = peaklist[:]  # 复制一个新列表,避免修改了原列表
    while i < len(newlist) - 1:
        j = i + 1
        while j < len(newlist):
            if between(newlist[i][0], newlist[j][0], size):
                average = format((newlist[i][0] + newlist[j][0]) / 2, '.4f')
                newlist[i] = (float(average), newlist[i][1] + newlist[j][1])
                newlist.remove(newlist[j])
                j = i
            j += 1
        i += 1
    return newlist


def addinformation(chem):
    '''传入一个STDARD_CHEM或者一个TEST_CHEM对象,为其添加W_list,distance,corr_distance'''
    chem.Weight = float(chem.Weight)
    if len(chem.Peak) == 0:
        return
    rpeak = combine_peak(chem.Peak)
    rpeak = sort_by_weight(rpeak)
    chem.relative_Peak = rpeak
    chem.W_list = W(rpeak)
    chem.distance = np.linalg.norm(np.array(chem.W_list))
    chem.corr_distance = np.linalg.norm(
        np.array(chem.W_list) - np.average(chem.W_list))


def cleandata(train_set):
    '''处理数据,将训练集分为已知和未知两部分,并为每个数据计算相对丰度、W_list和距离'''
    known = []
    unknown = []
    nopeak = []
    for chem in train_set:
        chem.next = []  # 添加一个参数，用于储存与其相连的后续节点
        if len(chem.Peak) == 0:
            nopeak.append(chem)
        elif chem.NAME == '':
            unknown.append(chem)
        else:
            known.append(chem)

    for kn in known:
        addinformation(kn)
    for un in unknown:
        un.marker = False  # False表示未知
        addinformation(un)
    return known, unknown
    # return known,unknown,nopeak


def abs_distance(standard, undertest):
    '''传入两个物质,返回两个物质的绝对值距离,注意该距离并非严格意义上的绝对值距离,且调整传入参数顺序将会得到不同的结果'''
    W1, W2 = insert_zero(standard, undertest)
    if len(W2) == 0:
        return 0
    absd = np.sum(np.abs(W1 - W2)) / np.sum(W2)
    return 1 / (1 + absd)


def euc_distance(standard, undertest):
    '''传入两个物质,返回这两个物质的欧式距离，注意该距离并非严格意义上的欧式距离,且调转参数传入顺序将得到不同结果'''
    W1, W2 = insert_zero(standard, undertest)
    if len(W2) == 0:
        return 0
    euc = np.sum((W1 - W2)**2) / np.sum(W2**2)
    return 1 / (1 + euc)


def correlation(testobj1, testobj2):
    '''传入两个物质,返回这两个物质的相关系数'''
    if len(testobj1.W_list) < 2 or len(testobj2.W_list) < 2:
        return 0
    W1, W2 = insert_zero(testobj1, testobj2)
    numer = np.sum((W1 - np.average(W1)) * (W2 - np.average(W2)))
    part1 = np.linalg.norm(
        np.array(testobj1.W_list) - np.average(testobj1.W_list))
    part2 = np.linalg.norm(
        np.array(testobj2.W_list) - np.average(testobj2.W_list))
    return abs(numer / (part1 * part2))


def Nei_similarity(testobj1, testobj2):
    '''改进的Nei系数法'''
    a1, a2 = insert_zero_absolute(testobj1, testobj2)
    cor = 2 / (len(testobj1.Peak) + len(testobj2.Peak))
    s = np.sum(np.abs((a1 - a2) / (a1 + a2)))
    return cor * (len(a1) - s)


def negative_index(standard, undertest):
    '''参数1为标准物质,参数2为待测物质;计算两个物质之间的新模型负指数,特别注意,如果更换两个参数的位置将会得到不一样的结果'''
    s, u = insert_zero_absolute(standard, undertest)
    if len(s) == 0:
        return 0
    sindex = np.sum(np.abs(s - u) / u)
    return np.exp(-1 * sindex / len(s))


def impro_similarity(standard, undertest):
    '''参数1为标准物质,参数2为待测物质;计算两个物质之间的改良程度相似度,特别注意,如果更换两个参数的位置将会得到不一样的结果'''
    s, u = insert_zero_absolute(standard, undertest)
    if len(s) == 0:
        return 0
    sindex = np.sum(np.abs(1 - s / u))
    return 1 - sindex / len(s)


def newimpro_similarity(standard, undertest):
    '''参数1为标准物质,参数2为待测物质;计算两个物质之间的新改良程度相似度,特别注意,如果更换两个参数的位置将会得到不一样的结果'''
    s, u = insert_zero_absolute(standard, undertest)
    if len(s) == 0:
        return 0
    sindex = np.sum((1 - s / u)**2)
    return 1 - np.sqrt(sindex / len(s))


def SS_similarity(standard, undertest):
    '''参数1为标准物质,参数2为待测物质;计算两个物质之间的Stein-Scott复合相似度,特别注意,如果更换两个参数的位置将会得到不一样的结果'''
    if len(standard.Peak) < 2:
        return 0
    n = 1 if standard.Peak[0][1] < standard.Peak[1][1] else -1
    s, u = insert_zero(standard, undertest)

    numer = np.sum(s * u)
    dot = 0 if numer == 0 else numer / (np.linalg.norm(s) * np.linalg.norm(u))
    # 计算点积,为了减少重复插0的不必要计算,这里没有直接调用dot_product()
    sr_part = (s[1:] * u[:-1]) / (s[:-1] * u[1:])
    sr = np.sum(sr_part**n)
    nx = len(undertest.Peak)
    nxy = len(u)
    return (nx * dot + nxy * sr) / (nxy + nx)


def dot_product(testobj1, testobj2):
    '''计算两个物质的点积,注意这里仅用到了丰度而没有使用到荷质比'''
    i1, i2 = insert_zero_absolute(testobj1, testobj2)
    numer = np.sum(i1 * i2)
    if numer == 0:
        return 0
    demon = np.linalg.norm(i1) * np.linalg.norm(i2)
    return numer / demon


def similarity_score(testobj1, testobj2):
    '''计算两个物质的相似性评分,注意仅用到了丰度而没有用到荷质比,该方法与点积类似,但对较低强度离子的权重更大'''
    i1, i2 = insert_zero_absolute(testobj1, testobj2)
    numer = np.sum(np.sqrt(i1 * i2))
    if numer == 0:
        return 0
    demon = np.sqrt(np.sum(i1) * np.sum(i2))
    return numer / demon


def similarity_index(testobj1, testobj2):
    '''计算两个物质的相似性因子,该算法使用到了丰度和荷质比两个因素'''
    W1, W2 = insert_zero(testobj1, testobj2)
    if len(W1) + len(W2) == 0:
        return 0
    numer = np.sum(np.abs(W1 - W2))
    demon = np.sum(W1 + W2)
    return 1 - numer / demon


def bonanza_similarity(testobj1, testobj2):
    '''计算两个物质的bonanza相似度'''
    one, two = 0, 0
    match1, match2, unmatch1, unmatch2 = [], [], [], []
    while one < len(testobj1.relative_Peak) and two < len(
            testobj2.relative_Peak):
        if between(testobj1.relative_Peak[one][0],
                   testobj2.relative_Peak[two][0]):
            w1 = testobj1.relative_Peak[one][0]**2 * np.sqrt(
                testobj1.relative_Peak[one][1])
            w2 = testobj2.relative_Peak[two][0]**2 * np.sqrt(
                testobj2.relative_Peak[two][1])
            match1.append(w1)
            match2.append(w2)
            one += 1
            two += 1
        elif testobj1.relative_Peak[one] < testobj2.relative_Peak[two]:
            w1 = testobj1.relative_Peak[one][0]**2 * np.sqrt(
                testobj1.relative_Peak[one][1])
            unmatch1.append(w1)
            one += 1
        else:
            w2 = testobj2.relative_Peak[two][0]**2 * np.sqrt(
                testobj2.relative_Peak[two][1])
            unmatch2.append(w2)
            two += 1
    if one < len(testobj1.relative_Peak):
        unmatch1.extend(
            [i[0]**2 * np.sqrt(i[1]) for i in testobj1.relative_Peak[one:]])
    if two < len(testobj2.relative_Peak):
        unmatch2.extend(
            [i[0]**2 * np.sqrt(i[1]) for i in testobj2.relative_Peak[two:]])
    mscore = np.sum(np.array(match1) * np.array(match2))
    uscore = np.sum(np.array(unmatch1)**2) + np.sum(np.array(unmatch2)**2)
    return mscore / (mscore + uscore)


def insert_zero(testobj1, testobj2):
    '''为W列表插入零向量,返回值为两个处理好的W列表'''
    one = 0
    two = 0
    copylist1 = testobj1.relative_Peak[:]
    copylist2 = testobj2.relative_Peak[:]
    W1 = testobj1.W_list[:]
    W2 = testobj2.W_list[:]
    while one < len(copylist1) and two < len(copylist2):
        if not between(copylist1[one][0], copylist2[two][0]):
            if copylist1[one][0] < copylist2[two][0]:
                del copylist1[one]
                del W1[one]
            else:
                del copylist2[two]
                del W2[two]
        else:
            one += 1
            two += 1
    if one < len(copylist1):
        W1 = W1[:one]
    if two < len(copylist2):
        W2 = W2[:two]
    return np.array(W1), np.array(W2)


def insert_zero_absolute(testobj1, testobj2):
    '''为物质的绝对丰度表插入零,返回两个处理好的峰度列表,仅包含对应的绝对含量而不带有荷质比'''
    one = 0
    two = 0
    copylist1 = testobj1.Peak[:]
    copylist2 = testobj2.Peak[:]
    while one < len(copylist1) and two < len(copylist2):
        if not between(copylist1[one][0], copylist2[two][0]):
            if copylist1[one][0] < copylist2[two][0]:
                del copylist1[one]
            else:
                del copylist2[two]
        else:
            one += 1
            two += 1
    if one < len(copylist1):
        copylist1 = copylist1[:one]
    if two < len(copylist2):
        copylist2 = copylist2[:two]
    copylist1 = np.array(copylist1).reshape(-1, 2)[:, 1]
    copylist2 = np.array(copylist2).reshape(-1, 2)[:, 1]
    return copylist1, copylist2


def dump_test_set():
    '''产生一个pickle文件,用于测试各种相似度公式的准确性,pickle文件是一个字典:
    {
        'standard':(一个已经被注释的TEST_CHEM对象,[标准谱库中与之同名的STD_CHEM对象列表]),
        'random':[未被注释的TEST_CHEM对象组成的列表,其峰的数量必大于2,对象随机抽取,列表长度与已注释对象数量相等]
    }'''
    from random import randint
    import pickle
    stdset = load_std_set()
    trainset = load_train_set()
    knownlist, unknownlist = [], []
    for t in trainset:
        if t.NAME != '':
            knownlist.append(t)
        elif len(t.Peak) > 2:
            unknownlist.append(t)
    standardlist = []
    for k in knownlist:
        thislist = [s for s in stdset if s.NAME == k.NAME]
        if len(thislist) == 0:
            continue
        standardlist.append((k, thislist))
    randomlist = []
    while len(randomlist) < len(standardlist):
        randomlist.append(unknownlist[randint(0, len(unknownlist) - 1)])
    pkldict = {'standard': standardlist, 'random': randomlist}
    with open("test_set.pkl", 'wb') as pklfile:
        pickle.dump(pkldict, pklfile)


def transform_fasttree(array, length):
    '''参数array为一个层次聚类好的列表,length为原始数据的个数,返回一个fasttree形式的字符串;
    注意,array中每行的参数为:[聚类簇1的序号, 聚类簇2的序号,聚类簇1 2之间的距离,聚类簇1的名字,聚类簇2的名字]
    如果参数对应顺序错误,将产生错误结果甚至程序暂停,对于新聚成的簇由于没有实际有意义的名字,可以随便写,这个不会被输出,但请保留这个参数的位置'''
    def recurse(index, array, length):
        nonlocal tree_string
        tree_string += '('
        left = array[index][0]
        right = array[index][1]
        if left > length:
            recurse(left - length - 1, array, length)
        if left <= length:
            tree_string += str(array[index][3])
        tree_string += ":{},".format(array[index][2])
        if right > length:
            recurse(right - length - 1, array, length)
        if right <= length:
            tree_string += str(array[index][4])
        tree_string += ":{})".format(array[index][2])

    tree_string = ''
    recurse(len(array) - 1, array, length)
    return tree_string + ';'


def mix_similarity(testobj1, testobj2):
    '''综合了cos,correlation,bonanza三种相似度的混合相似度,返回结果为三种相似度加和'''
    one, two = 0, 0
    match1, match2, unmatch1, unmatch2 = [], [], [], []
    while one < len(testobj1.relative_Peak) and two < len(
            testobj2.relative_Peak):
        if between(testobj1.relative_Peak[one][0],
                   testobj2.relative_Peak[two][0]):
            w1 = testobj1.relative_Peak[one][0]**2 * np.sqrt(
                testobj1.relative_Peak[one][1])
            w2 = testobj2.relative_Peak[two][0]**2 * np.sqrt(
                testobj2.relative_Peak[two][1])
            match1.append(w1)
            match2.append(w2)
            one += 1
            two += 1
        elif testobj1.relative_Peak[one] < testobj2.relative_Peak[two]:
            w1 = testobj1.relative_Peak[one][0]**2 * np.sqrt(
                testobj1.relative_Peak[one][1])
            unmatch1.append(w1)
            one += 1
        else:
            w2 = testobj2.relative_Peak[two][0]**2 * np.sqrt(
                testobj2.relative_Peak[two][1])
            unmatch2.append(w2)
            two += 1
    if one < len(testobj1.relative_Peak):
        unmatch1.extend(
            [i[0]**2 * np.sqrt(i[1]) for i in testobj1.relative_Peak[one:]])
    if two < len(testobj2.relative_Peak):
        unmatch2.extend(
            [i[0]**2 * np.sqrt(i[1]) for i in testobj2.relative_Peak[two:]])
    mscore = np.sum(np.array(match1) * np.array(match2))
    uscore = np.sum(np.array(unmatch1)**2) + np.sum(np.array(unmatch2)**2)
    match1, match2 = np.array(match1), np.array(match2)
    bonanza = mscore / (mscore + uscore)
    if testobj1.distance < 0.01 or testobj2.distance < 0.01:
        cos = 0
    else:
        cos = np.sum(match1 * match2) / (testobj1.distance * testobj2.distance)

    if len(match1) < 2 or len(match2) < 2:
        corr = 0
    else:
        corr_numer = np.sum(
            (match1 - np.average(match1)) * (match2 - np.average(match2)))
        corr = corr_numer / (testobj1.corr_distance * testobj2.corr_distance)
    return (bonanza + cos + corr) / 3
