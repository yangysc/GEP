# -*- coding:utf-8 -*-
from __future__ import division, print_function
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


# 操作符对应的操作数个数
# Q: sqrt,  I: if(if a =1, then b else c);   A: and,  O: or,  N : not
opt_arity = {'A': 2, 'O': 2, 'N': 1}   # function set
num_opt = len(opt_arity)  # number of function
terminal = {'a', 'b', 'c'}   # terminal
num_terminal = len(terminal)
T = list(range(len(terminal)))  # 实际染色体中变量取值
# symbols 中是 heads of the genes 可以取值的元素
symbols = list(opt_arity.keys())
symbols.extend(list(T))

# terminal 在染色体中为 0. 1. 2 ， 该下标 表示 保存在 下面的数组中的数值
symbol2values = [0] * len(T)
# the length of the head
h = 3
# maximum arity(单个函数最多的参数个数)
n_max = max(opt_arity.values())
# the length of the tail t is a function of h
t = h * (n_max - 1) + 1
# the length of the gene g
g = h + t
# the number of genes （一棵树等于 一个gene）
ngenes = 2
# 一条染色体长度
nvar = ngenes * g
# 存放操作符节点信息的队列
queue = []

labels = {}  # 图中的节点是字符的下标， labels记录下标的语义值(字符）

G = nx.DiGraph()  # 图

maxHeight = 99   # 允许的树的最大高度
# height_totaloffset = [0] * maxHeight

# def generateChromosome(x):
#     for i in range(h):
roots = [0] * ngenes

def translate(chromosome):
    # 每棵树的 root

    for i in range(ngenes):
        roots[i] = g * i
        parseOneGene(chromosome, roots[i])



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
        return symbol2values[int(labels[node])]
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
        elif labels[node] == 'Q':   # sqrt
            left = G[node].keys()[0]
            num = calculate(left)
            return np.sqrt(num)
        elif labels[node] == 'A':   # and
            left, right = G[node].keys()
            if calculate(left) == 1 and calculate(right) == 1:
                return 1
            else:
                return 0
        elif labels[node] == 'O':  # or
            left, right = G[node].keys()
            if calculate(left) == 1 or calculate(right) == 1:
                return 1
            else:
                return 0
        elif labels[node] == 'N':  # not
            left = G[node].keys()[0]
            if calculate(left) == 1:
                return 0
            else:
                return 1
        else:
            return 0


def init(x):
    for i in range(ngenes):
        for j in range(g):
            if j < num_opt:
                rnd_pos = np.random.randint(num_opt)
                x[i * g + j] = symbols[rnd_pos]
            else:
                rnd_pos = np.random.randint(num_terminal)
                x[i * g + j] = symbols[rnd_pos + num_opt]


def drawExprTree():
    g = pgv.AGraph()

    nodes = G.nodes()
    edges = G.edges()
    # print edges
    g.add_nodes_from(nodes)
    g.add_edges_from(edges)
    g.layout(prog="dot")

    for i in nodes:
        n = g.get_node(i)
        n.attr["label"] = str(labels[i]) + '(' + str(n) + ')'

    g.draw('tree.pdf')

def obj_eval(x):
    # truth table
    truth_table = np.array([[0, 0, 0, 0],
                            [0, 0, 1, 0],
                            [0, 1, 0, 0],
                            [0, 1, 1, 1],
                            [1, 0, 0, 0],
                            [1, 0, 1, 1],
                            [1, 1, 0, 1],
                            [1, 1, 1, 1]
                            ])
    f_val = 0
    translate(x)
    for rows in truth_table:
        symbol2values[:] = rows[:-1]
        value_exprTree = [0] * ngenes
        for i in range(ngenes):
            value_exprTree[i] = calculate(roots[i])

        finalOutcome = OR(value_exprTree[0], value_exprTree[1])
        # print(rows)
        # print('此时 每个基因结果分别 (%d, %d)' %(value_exprTree[0], value_exprTree[1]))
        # print('此时 y is %d' % finalOutcome)
        if finalOutcome == rows[-1]:
            f_val += 1
    return f_val


def OR(x1, x2):
    if x1 == 0 and x2 == 0:
        return 0
    else:
        return 1

if __name__ == '__main__':
    # c1 = "Q*-+2134"

    # translate(c1)
    # outcome = calculate(0)
    # expr = "np.sqrt((a-b) * (c + d))"
    # print "Algorithm gives %f" % outcome
    # trueVal = eval(expr)
    # print "True value is %f" % trueVal
    # x = ['a'] * nvar
    # symbol2values = []
    # init(x)
    # print(x)
    # translate(x)
    # drawExprTree()
    x1 = 'AAaccacNcaabab'
    x1 = x1.replace('a', '0')
    x1 = x1.replace('b', '1')
    x1 = x1.replace('c', '2')
    print(obj_eval(x1))

