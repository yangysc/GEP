# -*- coding:utf-8 -*-
# 回归浮点型参数
from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import pygraphviz as pgv
import time


# 遗传算法参数设置
maxIter = 300
popsize = 200

# 问题参数设置

all_best = []  # 存储每代最优个体


class QueueItem:
    def __init__(self, opt, idx, offset, uncle, nbrother, height):
        self.idx = idx  # 操作符所在下标
        self.offset = offset  # 节点与第一个孩子的相对距离， = 亲兄弟 + 堂兄弟 的个数）
        self.uncle = uncle  # 节点的叔叔们( 只记录操作符， terminal 没后代)
        self.nbrother = nbrother  # 右边的亲兄弟个数， 用于确定孩子有几个亲叔叔
        self.height = height     # 节点的高度


# 操作符对应的操作数个数
# Q: sqrt,  I: if(if a =1, then b else c);   A: and,  O: or,  N : not
opt_arity = {'+': 2, '-': 2, '*': 2, '/': 2}   # function set
num_opt = len(opt_arity)  # number of function
terminal = {'a', '?'}   # terminal
num_terminal = len(terminal)
# T = list(range(len(terminal)))  # 实际染色体中变量取值
# symbols 中是 heads of the genes 可以取值的元素
symbols = list(opt_arity.keys())
symbols.extend(list(terminal))
num_symbols = len(symbols)
# terminal 在染色体中为 0. 1. 2 ， 该下标 表示 保存在 下面的数组中的数值
symbol2values = [0] * num_terminal
# the length of the head
h = 7
# maximum arity(单个函数最多的参数个数)
n_max = max(opt_arity.values())
# the length of the tail t is a function of h
t = h * (n_max - 1) + 1
# the length of the gene g
g = h + t + t   # heads + tail + Dc
# the number of genes （一棵树等于 一个gene）
ngenes = 1
# 一条染色体长度
nvar = ngenes * g
# 存放操作符节点信息的队列
queue = []



# G = nx.DiGraph()  # 图

maxHeight = 99   # 允许的树的最大高度
# height_totaloffset = [0] * maxHeight

# def generateChromosome(x):
#     for i in range(h):
roots = [0] * ngenes

x_old = []
site_old = 0
len_old = 0
s_old = 0


# true function
def func(x):
    return x**2 / 2 + 3 * x  # 4 * x**4 + 3 * x**3 + 2 * x**2 + x


# 定义染色体
class Individual:
    def __init__(self):
        # 基因型
        self.x_var = ['0'] * nvar
        #  目标函数值（越小越好）
        self.y_obj = np.inf
        self.G = nx.DiGraph()  # 该染色体对应的表达树
        self.isValid = True
        self.constants = np.zeros(shape=(ngenes, 10), dtype=np.int)   # [[0] * 10 for i in range(ngenes)]    # 保存每个基因的浮点型参数
        self.numofconstant = [0] * ngenes   # 用于 obj_eval 时记录 每个基因已经访问了多少个 ?（参数）
        self.labels = {}  # 图中的节点是字符的下标， labels记录下标的语义值(字符）
        # init 个体基因型初始化及 估值
        self.init()

    # 个体初始化
    def init(self):

        for i in range(ngenes):
            # 初始化 基因型
            self.x_var[i * g + 0] = symbols[np.random.randint(num_opt)]
            for j in range(1, h):
                rnd_pos = np.random.randint(num_symbols)
                self.x_var[i * g + j] = symbols[rnd_pos]
            for j in range(h, h + t):
                rnd_pos = np.random.randint(num_terminal)
                self.x_var[i * g + j] = symbols[rnd_pos + num_opt]
            for j in range(h + t, g):  # Dc
                self.x_var[i * g + j] = str(np.random.randint(10))  # 注意用 char 型(与 c++ 保持一致， 虽然在python list中类型可以不一致）

            # 初始化 参数 constant (默认十个参数)
            for j in range(10):
                self.constants[i][j] = np.random.randint(10) # np.random.rand()

        self.obj_eval()

    # 适应度计算
    def obj_eval(self):
        self.G.clear()

        self.isValid = True
        self.labels = {}
        self.translate()

        # if self.isValid is False:
        #     self.y_obj = 0
        #     return
        # truth table
        truth_table = np.array([[6.9407, 44.909752],
                                [-7.8664, 7.3409245],
                                [-2.7861, -4.4771234],
                                [-5.0944, -2.3067443],
                                [9.4895, 73.493805],
                                [-9.6197, 17.410214],
                                [-9.4145, 16.072905],
                                [-0.1432, -0.41934688],
                                [0.9107, 3.1467872],
                                [2.1762, 8.8965232]
                                ])
        for ith, rows in enumerate(truth_table):
            truth_table[ith][0] = ith + 1
            truth_table[ith][-1] = func(rows[0])

        f_val = 0
        R = 100   # selection range
        for rows in truth_table:
            symbol2values[0] = rows[0]
            value_exprTree = [0] * ngenes
            for i in range(ngenes):

                value_exprTree[i] = self.calculate(roots[i])   # todo 这个parse 完生成树了之后，能不能递归一次就把所有的样本全计算了（数组计算）
                if self.isValid is False:
                    self.y_obj = 0
                    return

            finalOutcome = np.sum(value_exprTree)
            # print(rows)
            # print('此时 每个基因结果分别 ', value_exprTree)
            # print('此时 y is %f' % finalOutcome)

            f_val += R - np.abs(rows[-1] - finalOutcome)
        self.y_obj = f_val

    # 将染色体 p2 复制到自身
    def assignBy(self, p2):
        # 复制 基因型
        for i in range(nvar):
            self.x_var[i] = p2.x_var[i]
        # 复制 constants:
        for i in range(ngenes):
            for j in range(10):
                self.constants[i][j] = p2.constants[i][j]
        self.y_obj = p2.y_obj
        self.G.clear()

    # 输出在控制台
    def print(self):
        # x
        # print('x is:')
        # print(self.x_var)
        # y
        print('objective is %f' % self.y_obj)

    # 保存到本地
    def save(self):
        with open('outcome.txt', 'wr') as f:
            for i in range(nvar):
                print(i, sep=', ', file=f)
            print(self.y_obj, file=f)

    def translate(self):
        # 需要将字符序列变成 以下标记录的数组
        chromosome = list(range(len(self.x_var)))
        for node in chromosome:
            self.labels[node] = self.x_var[node]
        # 每棵树的 root
        self.numofconstant = [0] * ngenes  # 初始置零
        for i in range(ngenes):
            roots[i] = g * i
            self.parseOneGene(roots[i])


    def parseOneGene(self, root):

        # 记录本层中左边亲、堂兄弟中操作符 参数累计和
        height_totaloffset = [0] * maxHeight
        chromosome = list(range(len(self.x_var)))
        # 添加根节点
        pos = root
        gene = pos

        if self.labels[gene] not in opt_arity:  # 若根节点不是操作符
            if self.labels[gene] == "?":  # Dc

                #   取出 Dc
                currentGeneIdx = int(gene / g)
                # print(currentGeneIdx)
                DcIdx = self.x_var[currentGeneIdx * g + h + t + self.numofconstant[currentGeneIdx]]
                self.labels[gene] = DcIdx  # 直接改变'?' 的 labels
                self.numofconstant[currentGeneIdx] += 1
            self.G.add_node(gene)
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
            gene = pos  # gene = pos
            offset = item.offset  # 距第一个孩子的距离
            height = item.height  # 该节点在第几层
            uncle = item.uncle  # 该节点的叔叔们
            nbrother = item.nbrother  # 该节点的亲弟弟个数

            '''
            以下代码实现功能：
                            找到该节点的孩子，并将是操作符的孩子们添加到队列中

            '''
            max_arity = opt_arity[self.labels[gene]]  # 几个孩子
            child_idx = pos + offset + 1  # 第一个孩子的下标
            # 找出所有孩子
            children = chromosome[child_idx:  child_idx + max_arity]
            childIsOpt = False  # 记录孩子是否有操作符
            for child in children:
                if self.labels[child] in opt_arity:
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
                        nconsin += opt_arity[self.labels[u]]
                    cousin = chromosome[pos + 1 + nbrother : pos + 1 + + nbrother + nconsin]  # 不包括亲兄弟

                total_offset = 0  # 记录第一个孙子与该孩子的相对距离

                # 查找第一个孩子距第一个孙子的相对距离
                inc_ncousion = 0  # 由于自己【亲兄弟 + 堂兄弟】导致 产生的 自己孩子的堂兄弟的个数
                newuncle = []  # 记录自己的【亲兄弟 + 堂兄弟】中的【操作符】， 当作孩子们的叔叔们
                # 遍历自己亲兄弟
                for brother in chromosome[gene + 1: gene + 1 + nbrother]:
                    if self.labels[brother] in opt_arity:
                        newuncle.append(brother)
                        inc_ncousion += opt_arity[self.labels[brother]]
                # 若有堂兄弟
                if len(cousin) > 0:
                    # 遍历堂兄弟
                    # 找出堂兄弟中的操作符， 作为孩子们的叔叔
                    for item in cousin:
                        if self.labels[item] in opt_arity:
                            inc_ncousion += opt_arity[self.labels[item]]
                            newuncle.append(item)
                # 找自己孩子
                while current_arity < max_arity:
                    new_gene = chromosome[child_idx + current_arity]
                    # 如果孩子是操作符则入队
                    if self.labels[new_gene] in opt_arity:
                        newoffset = height_totaloffset[height + 2] + total_offset + (
                        max_arity - 1 - current_arity) + inc_ncousion
                        newitem = QueueItem(new_gene, child_idx + current_arity, newoffset, newuncle,
                                            (max_arity - 1 - current_arity), height + 1)
                        queue.append(newitem)
                        total_offset += opt_arity[self.labels[new_gene]]
                        # 添加边
                        self.G.add_edge(gene, new_gene)
                    else:  # 孩子是操作数
                        if self.labels[new_gene] == "?":  # Dc

                            #   取出 Dc
                            currentGeneIdx = int(new_gene / g)
                            # print(currentGeneIdx)
                            DcIdx = self.x_var[currentGeneIdx * g + h + t + self.numofconstant[currentGeneIdx]]
                            self.labels[new_gene] = DcIdx   # 直接改变'?' 的 labels
                            self.numofconstant[currentGeneIdx] += 1

                        self.G.add_edge(gene, new_gene)

                    current_arity += 1
                height_totaloffset[height + 2] += total_offset
            else:
                # 如果孩子中没有操作符，则不必计算孙子，只需要将孩子添加到树中
                for child in children:
                    if self.labels[child] == "?":  # Dc

                        #   取出 Dc
                        currentGeneIdx = int(child / g)
                        # print(currentGeneIdx)
                        DcIdx = self.x_var[currentGeneIdx * g + h + t + self.numofconstant[currentGeneIdx]]
                        self.labels[child] = DcIdx  # 直接改变'?' 的 labels
                        self.numofconstant[currentGeneIdx] += 1
                    self.G.add_edge(gene, child)
            # 队列中没有元素(操作符), 说明最后一层均为 操作数， 退出循环
            if len(queue) < 1:
                return

                # 递归计算树的算术值

    def calculate(self, node):
        if self.labels[node] not in opt_arity:
            if '0' <= self.labels[node] <= '9':
                # 该节点位于第几个基因
                currentGeneIdx = int(node / g)
                # print(currentGeneIdx)
                DcIdx = int(self.labels[node])
                c = self.constants[currentGeneIdx][DcIdx]
                return c
            else:
                try:
                    return float(symbol2values[ord(self.labels[node]) - ord('a')])
                except IndexError as ie:
                    print(ie.message)
                    print("node is %d, label is %s " %(node, self.labels[node]))
                    print(self.labels)
                    print(" (%d, %d)" % (len(self.x_var), nvar))
                    exit(-10)
        else:
            if self.labels[node] == '+':
                left, right = self.G[node].keys()
                if left > right:
                    temp = right
                    right = left
                    left = temp
                leftValue = self.calculate(left)
                rightValue = self.calculate(right)
                return leftValue + rightValue
            elif self.labels[node] == '-':
                left, right = self.G[node].keys()
                if left > right:
                    temp = right
                    right = left
                    left = temp
                leftValue = self.calculate(left)
                rightValue = self.calculate(right)
                return leftValue - rightValue

            elif self.labels[node] == '*':
                left, right = self.G[node].keys()
                if left > right:
                     temp = right
                     right = left
                     left = temp
                leftValue = self.calculate(left)
                rightValue = self.calculate(right)
                return leftValue * rightValue

            elif self.labels[node] == '/':
                left, right = self.G[node].keys()

                # .keys返回的两个孩子不分左右，但是程序里左孩子的 index 均小于右孩子
                if left > right:
                    temp = right
                    right = left
                    left = temp
                # 被零除返回 1
                rightValue = self.calculate(right)
                if rightValue == 0:
                    self.isValid = False
                    return 1
                else:
                    return self.calculate(left) / rightValue
            elif self.labels[node] == 'Q':  # sqrt
                left = self.G[node].keys()[0]

                num = self.calculate(left)
                return np.sqrt(num)
            elif self.labels[node] == 'A':  # and
                left, right = self.G[node].keys()
                if left > right:
                    temp = right
                    right = left
                    left = temp
                if self.calculate(left) == 1 and self.calculate(right) == 1:
                    return 1
                else:
                    return 0
            elif self.labels[node] == 'O':  # or
                left, right = self.G[node].keys()
                if left > right:
                    temp = right
                    right = left
                    left = temp
                if self.calculate(left) == 1 or self.calculate(right) == 1:
                    return 1
                else:
                    return 0
            elif self.labels[node] == 'N':  # not
                left = self.G[node].keys()[0]

                if self.calculate(left) == 1:
                    return 0
                else:
                    return 1
            else:
                return 0

    def drawExprTree(self, src):
        g = pgv.AGraph()

        nodes = self.G.nodes()
        edges = self.G.edges()
        # print edges
        g.add_nodes_from(nodes)
        g.add_edges_from(edges)
        g.layout(prog="dot")

        for i in nodes:
            n = g.get_node(i)
            n.attr["label"] = str(self.labels[i]) + '(' + str(n) + ')'
        g.draw(src + '.pdf')


# GEP 算法程序
class GEP:
    def __init__(self):
        # 产生种群
        self.population = [Individual() for i in range(popsize)]
        self.offspring = [Individual() for i in range(2 * popsize)]
        self.bestIndividual = Individual()
        self.iter = 0  # 当前是第几次循环， 用于 NUM

        self.CXPB = 0.7  # 交叉概率
        self.MUTPB = 2 / nvar  # 变异概率 (equivalent to two one-point mutations per chromosome
        self.InvPB = 0.1   # inversion rate
        self.ISPB = 0.1   # IS transposition rate
        self.RISPB = 0.1   # root transposition
        self.GENEPB = 0.1  # gene transposition rate

    # 二进制随机联赛(用于从父代中选择子代）
    def bin_tournament(self):
        p1 = np.random.randint(popsize)
        p2 = np.random.randint(popsize)
        while p2 == p1:
            p2 = np.random.randint(popsize)
        if self.population[p1].y_obj > self.population[p2].y_obj:
            return p1
        else:
            return p2

    # 选择
    def select(self):
        # 锦标赛选择(从父代+子代中选出新的父代)
        for i in range(popsize):
            p1 = np.random.randint(2 * popsize)
            p2 = np.random.randint(2 * popsize)
            while p2 == p1:
                p2 = np.random.randint(2 * popsize)
            if self.offspring[p1].y_obj > self.offspring[p2].y_obj:
                self.population[i].assignBy(self.offspring[p1])
            else:
                self.population[i].assignBy(self.offspring[p2])

        # 精英保留策略，用上一代的最优的替代本次种群的最差的( 目标函数越小越好)
        minV = self.population[0].y_obj
        min_idx = 0
        for i in range(1, popsize):
            if minV < self.population[i].y_obj:
                minV = self.population[i].y_obj
                min_idx = i
        self.population[min_idx].assignBy(self.bestIndividual)

    # 找出 Population 中 最优染色体（y_obj 最大）, 保存之
    def storeBest(self):
        best_y = -np.inf
        best_idx = -1
        for i in range(popsize):  # 居然把self.popsize 写成了self.nvar , 打你哭 （感谢我机智的debug ability）
            if best_y < self.population[i].y_obj:
                best_y = self.population[i].y_obj
                best_idx = i
        self.bestIndividual.assignBy(self.population[best_idx])

    # Mutation
    def mutate(self, x):
        # 变异位数 ( 两位 )
        numofMutation = int(self.MUTPB * nvar)
        for i in range(numofMutation):
            # 先选择哪条基因
            rnd_gene = np.random.randint(ngenes)
            # 选择基因中位置 (heads + tail)
            rnd_pos = np.random.randint(h + t)
            if rnd_pos < h:   # head
                x[rnd_gene * g + rnd_pos] = np.random.choice(symbols)
            else:  # tail
                rnd_symbol = np.random.randint(num_terminal)    # 随机从terminal 中选择一个下标
                x[rnd_gene * g + rnd_pos] = symbols[rnd_symbol + num_opt]


    # Mutation of Dc
    def mutateDc(self, x):
        # 变异位数 ( 两位 )
        numofMutation = int(self.MUTPB * nvar)
        for i in range(numofMutation):
            # 先选择哪条基因
            rnd_gene = np.random.randint(ngenes)
            # 选择基因中位置 (heads + tail)
            rnd_pos = np.random.randint(t)

            x[rnd_gene * g + h + t + rnd_pos] = str(np.random.randint(10))
    # Inversion
    def inversion(self, x):
        # select a gene
        gene = np.random.randint(ngenes)
        # the start and termination points of the sequence
        startP = np.random.randint(h - 1)
        endP = np.random.randint(startP, h)
        midP = int((startP + endP) / 2)
        # 如果选择的区间长度为奇数， 如[0, 3], 说明有偶数个元素
        if (endP - startP) % 2 == 1:
            midP += 1


        # reverse
        for i in range(0, midP - startP):
            temp = x[gene * g + startP + i]
            x[gene * g + startP + i] = x[gene * g + endP - i]
            x[gene * g + endP - i] = temp

    def inversion_Dc(self, x):  # inversion Dc
        # select a gene
        gene = np.random.randint(ngenes)
        # the start and termination points of the sequence
        startP = np.random.randint(t - 1)  # len(Dc) = t
        endP = np.random.randint(startP, t)

        # test
        # gene = 0
        # startP = 0
        # endP = 3

        midP = int((startP + endP) / 2)
        # 如果选择的区间长度为奇数， 如[0, 3], 说明有偶数个元素
        if (endP - startP) % 2 == 1:
            midP += 1

        # reverse
        for i in range(0, midP - startP):
            temp = x[gene * g + startP + h + t + i]
            x[gene * g + startP + h + t + i] = x[gene * g + + h + t + endP - i]
            x[gene * g + + h + t + endP - i] = temp

    #  transposition of IS elements
    def transpositionIS(self, x):

        # select a gene
        gene = np.random.randint(ngenes)


        # print(len)
        # insert site ( must in heads)
        site = np.random.randint(1, h - 1)      # 起点必须 > 0 !!!!
        # the start point of the sequence
        startP = np.random.randint(site, h + t - 1)  # 起点必须 > 0 !!!!
        # the length of the sequence
        length = np.random.randint(1, min(h - site, h + t - startP))  # length + site 不能超过 h...and length + startP <= g
        '''debug'''

        # x_old[:] = x[:]
        # global len_old
        # global site_old
        # global s_old
        # site_old = site
        # len_old = length
        # s_old = startP

        ### 测试
        # gene = 0
        # startP = 14
        # len = 3
        # site = 1
        # 先将 site + len: h 的元素向后移动 len 个单位
        backup = x[gene * g + startP: gene * g + startP + length][:]
        for i in range(h - length - site):
            x[gene * g + site + length + i] = x[gene * g + site + i]
        # print('done')
        # 复制 transposition 的基因
        x[gene * g + site: gene * g + site + length] = backup

    #  transposition of RIS elements
    def transpositionRIS(self, x):

        # select a gene
        gene = np.random.randint(ngenes)

        # insert site ( must in heads)
        startP = np.random.randint(1, h)  # 随机选择一个点
        # scan the gene until a function is found
        while x[startP] not in opt_arity and startP < h + t:
            startP += 1
        if startP >= h + t:
            return



        # the length of the sequence
        length = np.random.randint(1, max(h - startP, 2))  # length + site 不能超过 h...and length + startP <= g

        x_old[:] = x[:]
        global len_old
        global site_old
        global s_old
        site_old = 0
        len_old = length
        s_old = startP

        # 先将 site + len: h 的元素向后移动 len 个单位
        backup = x[gene * g + startP: gene * g + startP + length][:]
        for i in range(h - length):
            x[gene * g + length + i] = x[gene * g + i]
        # print('done')
        # 复制 transposition 的基因
        x[gene * g: gene * g + length] = backup


    def transpositionGene(self, x):
        # select a gene
        gene = np.random.randint(ngenes)
        # transpose it to the begining of the chromosome
        backup = x[gene * g: gene * g + g][:]
        for i in range(g, gene * g + g):
            x[i] = x[i - g]
        x[: g] = backup

        #  transposition of IS elements

    def transpositionDc(self, x):

        # select a gene
        gene = np.random.randint(ngenes)

        # print(len)
        # insert site ( must in heads)
        site = np.random.randint(0, t - 2)  # 起点必须 > 0 !!!!
        # the start point of the sequence
        startP = np.random.randint(site, t - 1)  # 起点必须 > 0 !!!!
        # the length of the sequence
        length = np.random.randint(1, min(2, t - startP))  # length + site 不能超过 h...and length + startP <= g
        '''debug'''

        # x_old[:] = x[:]
        # global len_old
        # global site_old
        # global s_old
        # site_old = site
        # len_old = length
        # s_old = startP

        ### 测试
        # gene = 1
        # startP = 3
        # length = 3
        # site = 0

        # 先将 site + len: h 的元素向后移动 len 个单位
        backup = x[gene * g + h + t + startP: gene * g + startP + h + t + length][:]
        for i in range(t - length - site):
            x[gene * g + h + t + site + length + i] = x[gene * g + h + t + site + i]
        # print('done')
        # 复制 transposition 的基因
        x[gene * g + h + t + site: gene * g + h + t + site + length] = backup


    def crossover(self, p1, p2, child):
        # # 单点交叉
        # # 随机选择交叉点
        # rnd_pos = np.random.randint(1, nvar - 1)
        # for i in range(rnd_pos):
        #     child[i] = p1[i]
        # for i in range(rnd_pos, nvar):
        #     child[i] = p2[i]

        # 两点交叉
        # 随机选择两个交叉点
        pos1 = np.random.randint(1, nvar - 1)
        pos2 = np.random.randint(pos1, nvar - 1)
        for i in range(pos1):
            child[i] = p1[i]
        for i in range(pos1, pos2):
            child[i] = p2[i]
        for i in range(pos2, nvar):
            child[i] = p1[i]

    def evolution(self):
        # generate offsprings
        for i in range(popsize):
            # 找出用来和 i 重组的个体
            p1 = self.bin_tournament()
            # todo  p2随机产生，还是锦标赛产生？

            self.offspring[i].assignBy(self.population[p1])

            # # 锦标赛产生
            p2 = self.bin_tournament()
            while p2 == p1:
                p2 = self.bin_tournament()
            if np.random.rand() < self.CXPB:
                self.crossover(self.population[p1].x_var, self.population[p2].x_var, self.offspring[i].x_var)
            else:
                self.offspring[i].assignBy(self.population[p1])
            self.mutate(self.offspring[i].x_var)
            self.mutateDc(self.offspring[i].x_var)

            #
            # transpositi
            if np.random.rand() < self.InvPB:
                self.inversion(self.offspring[i].x_var)
            if np.random.rand() < self.ISPB:
                self.inversion_Dc(self.offspring[i].x_var)
            if np.random.rand() < self.ISPB:
                self.transpositionIS(self.offspring[i].x_var)
            if np.random.rand() < self.RISPB:
                self.transpositionRIS(self.offspring[i].x_var)

            if np.random.rand() < self.GENEPB:
                self.transpositionGene(self.offspring[i].x_var)

            if np.random.rand() < 0.3: # self.ISPB:
                self.transpositionDc(self.offspring[i].x_var)


            # try:
            self.offspring[i].obj_eval()
            # except ValueError as VE:
            #     print(VE.message)
            #     print(self.offspring[i].x_var)
            #     print("Old")
            #     print(x_old)
            #     print("Trans site is %d, len is %d, start point is %d" % (site_old, len_old, s_old))
            #     exit(-10)

            # 将父代添加至后代中
            self.offspring[i + popsize].assignBy(self.population[i])
        # select next generation from parental population and offsprings

        self.select()


    # 迭代主程序
    def run(self):
        # # 找出当前最好的个体
        self.storeBest()
        all_best.append(self.bestIndividual.y_obj)
        # evolution
        for i in range(maxIter):
            self.iter = i
            self.evolution()
            self.storeBest()
            all_best.append(self.bestIndividual.y_obj)

            print('Running the %d iteration!!' % i)
            print(self.bestIndividual.x_var)
            self.bestIndividual.print()

        self.bestIndividual.print()
        print(self.bestIndividual.x_var)
        print("constants are:")
        print("%s" % ", ".join([str(item) for item in self.bestIndividual.constants[0]]))
        self.bestIndividual.obj_eval()
        self.bestIndividual.print()
        self.bestIndividual.drawExprTree('bestIndividual')

#
# def drawExprTree():
#     g = pgv.AGraph()
#
#     nodes = G.nodes()
#     edges = G.edges()
#     # print edges
#     g.add_nodes_from(nodes)
#     g.add_edges_from(edges)
#     g.layout(prog="dot")
#
#     for i in nodes:
#         n = g.get_node(i)
#         n.attr["label"] = str(labels[i]) + '(' + str(n) + ')'
#
#     g.draw('tree.pdf')
#
# def obj_eval(x):
#     # truth table
#     truth_table = np.array([[0, 0, 0, 0],
#                             [0, 0, 1, 0],
#                             [0, 1, 0, 0],
#                             [0, 1, 1, 1],
#                             [1, 0, 0, 0],
#                             [1, 0, 1, 1],
#                             [1, 1, 0, 1],
#                             [1, 1, 1, 1]
#                             ])
#     f_val = 0
#     translate(x)
#     for rows in truth_table:
#         symbol2values[:] = rows[:-1]
#         value_exprTree = [0] * ngenes
#         for i in range(ngenes):
#             value_exprTree[i] = calculate(roots[i])
#
#         finalOutcome = OR(value_exprTree[0], value_exprTree[1])
#         # print(rows)
#         # print('此时 每个基因结果分别 (%d, %d)' %(value_exprTree[0], value_exprTree[1]))
#         # print('此时 y is %d' % finalOutcome)
#         if finalOutcome == rows[-1]:
#             f_val += 1
#     return f_val


def OR(x):
    if np.any(np.array(x) == 1):
        return 1
    else:
        return 0


def main():
    np.random.seed()
    # new instance of RCGA
    gep = GEP()
    gep.run()
    plt.plot(all_best)
    plt.xlabel('iteration')
    plt.ylabel('objective')
    plt.grid()
    plt.show()

if __name__ == '__main__':
    main()
    # x1 = "*+/+/a?a??a407309+++*+??a???737256"
    # x1 = ['/', 'a', '/', '?', '+', '?', '-', '?', 'a', 'a', 'a', 'a', 'a', 'a', '?', '9', '0', '3', '6', '4', '9', '3', '6']
    # x1 = "*+/+a++??aa????20925007"
    # # x1 = x1.replace('a', '0')
    # # x1 = x1.replace('b', '1')
    # # x1 = x1.replace('c', '2')
    # ind = Individual()
    # ind.x_var[:] = list(x1)
    # ind.constants[0][:] = [6, 0, 0, 7, 0, 7, 0, 5, 8, 4]
    # print(ind.x_var)
    # ind.obj_eval()
    # print(ind.y_obj)
    # ind.drawExprTree('before')
    # gep = GEP()
    # gep.inversion_Dc(ind.x_var)
    # ind.obj_eval()
    # print(ind.x_var)
    # print(ind.y_obj)
    # ind.drawExprTree('after')

