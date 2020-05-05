class TreeNode(object):
    max_length = 0

    def __init__(self):
        self.father = None
        self.left = None
        self.right = None
        self.name = ''
        self.value = 0
        self.direction = 'left'


def transform_binary(nwk_string):
    '''将传入的字符串nwk文件转为二叉树'''
    node_pointer = TreeNode()
    root_node = node_pointer
    TreeNode.max_length = find_max_length(nwk_string)
    i = 1  # 由于前面已经建立了一个根节点,因此从1开始
    nodelist = [root_node]  # 节点列表,列表中存储了所有节点的地址
    meet_confidence = False
    while i < len(nwk_string):
        char = nwk_string[i]
        i += 1
        if char == '(':
            node = TreeNode()
            nodelist.append(node)
            if node_pointer.direction == 'left':
                node_pointer.left = node
            elif node_pointer.direction == 'right':
                node_pointer.right = node
            node.father = node_pointer
            node_pointer = node
        elif char == ':':
            nextchar = nwk_string[i]
            while ord('0') <= ord(nextchar) <= ord('9') or nextchar == '.':
                char += nextchar
                i += 1
                nextchar = nwk_string[i]
            node_pointer.value = float(char[1:])
            node_pointer = node_pointer.father
        elif char == ',':
            if node_pointer.direction == 'right':
                raise AssertionError  # 如果出现这个错误说明这个树不是二叉的,或者nwk格式文件有误
            node_pointer.direction = 'right'
        elif char == ')':
            meet_confidence = True
        elif char == ';':
            return nodelist  # 返回在这里
        elif char == ' ' or char == '\n' or char == '\t' or char == "'":
            continue
        elif meet_confidence:
            nextchar = nwk_string[i]
            while ord('0') <= ord(nextchar) <= ord('9') or nextchar == '.':
                char += nextchar
                i += 1
                nextchar = nwk_string[i]
            meet_confidence = False
        else:
            node = TreeNode()
            nodelist.append(node)
            if node_pointer.direction == 'left':
                node_pointer.left = node
            elif node_pointer.direction == 'right':
                node_pointer.right = node
            node.father = node_pointer
            node_pointer = node
            while nwk_string[i] != ':':
                nextchar = nwk_string[i]
                char += nextchar
                i += 1
            node.name = char


def find_max_length(nwk_string):
    left_branket, right_branket, max_length = 0, 0, 0
    for s in nwk_string:
        if s == '(':
            left_branket += 1
            max_length = max(max_length, left_branket - right_branket)
        elif s == ')':
            right_branket += 1
    return max_length


def find_all_children(root_node):
    '''传入一个父节点,该方法能返回该父节点下所有非空子节点的值'''
    def traverse(node):
        nonlocal result_list
        if node != None:
            traverse(node.left)
            if node.name != '':
                result_list.append(node.name)
            traverse(node.right)

    result_list = []
    traverse(root_node)
    return result_list


def all_in(standard, test):
    '''检测standard中所有的元素是否都在test中存在'''
    for s in standard:
        if s not in test:
            return False
    return True


def compare_timetree(namelist, cladelist, binarylist):
    '''该函数是对TimeTree结果进行打分的'''
    result_list, new_namelist, new_cladelist = [], [], []
    for i in range(len(namelist)):
        if namelist[i] not in new_namelist:
            new_namelist.append(namelist[i])
            new_cladelist.append(cladelist[i])  # 去重
    namelist, cladelist = new_namelist, new_cladelist
    for i in binarylist:
        if i.name != '':
            i.name = i.name.replace("_", " ")
    treename = [i.name for i in binarylist if i.name != '']
    index = 0
    while index < len(namelist):
        if namelist[index] not in treename:
            del namelist[index]
            del cladelist[index]
        else:
            index += 1
    for clade in set(cladelist):
        thislist = [
            namelist[i] for i in range(len(cladelist)) if cladelist[i] == clade
        ]
        for node in binarylist:
            if node.name == thislist[0]:  #首先找到第一个节点
                break
        children = find_all_children(node)
        while not all_in(thislist, children):
            node = node.father
            children = find_all_children(node)
        minscore = len(thislist) / len(namelist)  # 假设该clade分散在整个树中,也能得到这个分数
        score = (len(thislist) / len(children) - minscore) / (1 - minscore)
        result_list.append((clade, score))
    return result_list


def compare_mytree(binarylist):
    '''该函数是对自己建树的结果进行打分的'''
    result_list = []
    namelist = [i.name for i in binarylist if i.name != ""]
    cladelist = list(set(namelist))
    for clade in cladelist:
        totalnum = 0  # 统计这个clade共有多少个
        for n in namelist:
            if n == clade:
                totalnum += 1
        for node in binarylist:
            if node.name == clade:
                break  # 先找到第一个node
        while True:
            thisnum = 0
            children = find_all_children(node)
            for c in children:
                if c == clade:
                    thisnum += 1
            if thisnum >= totalnum:
                break
            node = node.father
            if node.father is None:
                a = 1
        minscore = totalnum / len(namelist)
        score = (totalnum / len(children) - minscore) / (1 - minscore)
        result_list.append((clade, score))
    return result_list


def main():
    namelist, cladelist = [], []
    with open("data/Name and Clade.txt", 'rt') as txtfile:
        # 该文件中存储了标记好的物种与其对应的clade
        for lines in txtfile:
            lines = lines.replace("\n", '')
            lines = lines.split("\t")
            namelist.append(lines[0])
            cladelist.append(lines[1])

    with open("data/timetree.nwk", 'rt') as nwkfile:
        # 该文件为TimeTree建树结果
        timetree_list = transform_binary(nwkfile.read())
    result_list = compare_timetree(namelist, cladelist, timetree_list)
    with open("data/timetree_score.txt", 'wt') as txtfile:
        for r in result_list:
            txtfile.write("{}\t{}\n".format(r[0], r[1]))
    # 以上部分为对TimeTree结果进行打分

    with open("600sample/clade_result.nwk", 'rt') as nwkfile:
        # 该文件为construct_phylogenic_tree.py中输出的结果
        mytree_list = transform_binary(nwkfile.read())
    result_list = compare_mytree(mytree_list)
    with open("600sample/mytree_score.txt", 'wt') as txtfile:
        for r in result_list:
            txtfile.write("{}\t{}\n".format(r[0], r[1]))


if __name__ == '__main__':
    main()