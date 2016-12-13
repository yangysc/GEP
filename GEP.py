# -*- coding:utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import pygraphviz as pgv
import time


class QueueItem:
    def __init__(self, opt, idx, offset, uncle, nbrother, height):
        self.idx = idx  # 操作符所在下标
        self.offset = offset  # 节点与第一个孩子的相对距离， = 亲兄弟 + 堂兄弟 的个数）
        self.uncle = uncle  # 节点的叔叔们( 只记录操作符， terminal 没后代)
        self.nbrother = nbrother  # 右边的亲兄弟个数， 用于确定孩子有几个亲叔叔
        self.height = height     # 节点的高度

# the length of the head
h = 5

# 操作符对应的操作数个数
# Q: sqrt,  I: if(if a =1, then b else c);   A: and,  O: or,  N : not
opt_arity = {'Q': 1, '+': 2, '-': 2, '/': 2, '*': 2,   'I': 3, 'A': 2, 'O': 2, 'N': 1}
# maximum arity(单个函数最多的参数个数)
n_max = max(opt_arity.values())
# the length of the tail t is a function of h
t = h * (n_max - 1) + 1
# the length of the gene g
g = h + t
# the number of genes （一棵树等于 一个gene）
ngenes = 3
# 存放操作符节点信息的队列
queue = []

labels = {}  # 图中的节点是字符的下标， labels记录下标的语义值(字符）

G = nx.DiGraph()  # 图

maxHeight = 99   # 允许的树的最大高度
# height_totaloffset = [0] * maxHeight

# def generateChromosome(x):
#     for i in range(h):


def translate(chromosome):
    # 每棵树的 root
    # roots = [0] * ngenes
    for i in range(ngenes):
        # roots[i] = g * i
        parseOneGene(chromosome, g * i)
        # height_totaloffset = [0] * maxHeight


def parseOneGene(chromo, root):

    # 需要将字符序列变成 以下标记录的数组
    chromosome = list(range(len(chromo)))
    for node in chromosome:
        labels[node] = chromo[node]
    # 记录本层中左边亲、堂兄弟中操作符 参数累计和
    height_totaloffset = [0] * maxHeight
    # 添加根节点
    pos = root
    gene = chromosome[pos]

    if labels[gene] not in opt_arity:  # 若根节点不是操作符
        G.add_node(gene)
        return
    else:  # 是操作符则将根节点添加到队列中
        offset = 0
        height = 0
        uncle = []
        nbro = 0
        queue.append(QueueItem(gene, pos, offset, uncle, nbro, height))
    # 循环添加其他节点
    while 1:
        # 弹出节点（该节点(队列中)一定是操作符）
        item = queue.pop(0)
        pos = item.idx  # 下标
        gene = pos      # gene = pos
        offset = item.offset  # 距第一个孩子的距离
        height = item.height  # 该节点在第几层
        uncle = item.uncle    # 该节点的叔叔们
        nbrother = item.nbrother    # 该节点的亲弟弟个数

        '''
        以下代码实现功能：
                        找到该节点的孩子，并将是操作符的孩子们添加到队列中

        '''
        max_arity = opt_arity[labels[gene]]  # 几个孩子
        child_idx = pos + offset + 1  # 第一个孩子的下标
        # 找出所有孩子
        children = chromosome[child_idx:  child_idx + max_arity ]
        childIsOpt = False  # 记录孩子是否有操作符
        for child in children:
            if labels[child] in opt_arity:
                childIsOpt = True
                break
        # 如果孩子中有操作符 则继续以下操作
        if childIsOpt:
            current_arity = 0

            '''
                找到节点的堂兄弟，用于确定本层 右侧的元素（操作符 + 操作数)个数

            '''
            cousin = []
            nconsin = 0
            # 若有叔叔，则遍历叔叔们，找到 堂兄弟个数(不包括自己的亲兄弟, 堂兄弟中也包括操作数)
            if len(uncle) > 0:
                for u in uncle:
                    nconsin += opt_arity[labels[u]]
                cousin = chromosome[pos + max_arity: pos + max_arity + nconsin]  # 不包括亲兄弟

            total_offset = 0               # 记录第一个孙子与该孩子的相对距离

            # 查找第一个孩子距第一个孙子的相对距离
            inc_ncousion = 0  # 由于自己【亲兄弟 + 堂兄弟】导致 产生的 自己孩子的堂兄弟的个数
            newuncle = []     # 记录自己的【亲兄弟 + 堂兄弟】中的【操作符】， 当作孩子们的叔叔们
            # 遍历自己亲兄弟
            for brother in chromosome[gene + 1: gene + 1 + nbrother]:
                if labels[brother] in opt_arity:
                    newuncle.append(brother)
                    inc_ncousion += opt_arity[labels[brother]]
            # 若没有堂兄弟
            if len(cousin) > 0:
                pass
            else:   # 遍历堂兄弟
                # 找出堂兄弟中的操作符， 作为孩子们的叔叔
                tempuncle = [item for item in cousin if labels[item] in opt_arity]
                newuncle.extend(tempuncle)
            # 找自己孩子
            while current_arity < max_arity:
                new_gene = chromosome[child_idx + current_arity]
                # 如果孩子是操作符则入队
                if labels[new_gene] in opt_arity:
                    newoffset = height_totaloffset[height + 1] + total_offset + (max_arity - 1 - current_arity) + nconsin + inc_ncousion
                    newitem = QueueItem(new_gene, child_idx + current_arity, newoffset,  newuncle, (max_arity - 1 - current_arity), height + 1)
                    queue.append(newitem)
                    total_offset += opt_arity[labels[new_gene]]
                # 添加边
                G.add_edge(gene, new_gene)
                current_arity += 1
            height_totaloffset[height + 1] += total_offset
        else:
            # 如果孩子中没有操作符，则不必计算孙子，只需要将孩子添加到树中
            for child in children:
                G.add_edge(gene, child)
        # 队列中没有元素(操作符), 说明最后一层均为 操作数， 退出循环
        if len(queue) < 1:
            return G


# 递归计算树的算术值
def calculate(node):
    if labels[node] not in opt_arity:
        return float(labels[node])
    else:
        if labels[node] == '+':
            left, right = G[node].keys()
            return calculate(left) + calculate(right)
        elif labels[node] == '-':
            left, right = G[node].keys()
            return calculate(left) - calculate(right)
        elif labels[node] == '*':
            left, right = G[node].keys()
            return calculate(left) * calculate(right)
        elif labels[node] == '/ ':
            left, right = G[node].keys()
            # 被零除返回 1
            if right == 0:
                return 1
            else:
                return calculate(left) / calculate(right)
        elif labels[node] == 'Q':
            left = G[node].keys()[0]
            num = calculate(left)
            return np.sqrt(num)
        else:
            return 0


if __name__ == '__main__':
    c1 = "Q*-+2134"
    c2 = "Q*b**+baQba"

    c4 = '+++++++00000000'
    c5 = '*b+a-+Qab+//+b+babbabbbababbaaa'
    c6 = '++++++bababcd'
    c7 = 'IaIcaIcabc'
    c8 = 'NIAbObbaaaabaabb'
    c9 = 'AOaabaaaabNabaaaaaabINNbababaa'
    c10 = '*Qb+*/bbbabab-a+QbQbbababa/ba-/*bbaaaaa'
    c11 = '-/dac/dacaccd//-aacbbbabcd-d/+c*def'
    c12 = 'QaQ+-Qbbaaaba+Q+ab+abababa*-**b+aabbaba'
    c13 = 'IOaIAcbaaacaacacAOcaIccabcbccbacIONAAbbbbacbcbbc'
    translate(c13)
    # outcome = calculate(0)
    # expr = "np.sqrt((a-b) * (c + d))"
    # print "Algorithm gives %f" % outcome
    # trueVal = eval(expr)
    # print "True value is %f" % trueVal
    print nx.info(G)
    g = pgv.AGraph()

    nodes = G.nodes()
    edges = G.edges()
    # print edges
    g.add_nodes_from(nodes)
    g.add_edges_from(edges)
    g.layout(prog="dot")

    for i in nodes:
        n = g.get_node(i)
        n.attr["label"] = labels[i] + '('+ str(n) + ')'

    nodes = G.nodes()
    g.draw('tree6.pdf')
