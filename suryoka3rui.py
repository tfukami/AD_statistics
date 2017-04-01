import numpy as np
from matplotlib.font_manager import FronProperties
from matplotlib.pyplot import figure, show
import matplotlib.pyplot as plt
import csv

csvfile = open('data.txt')
dlist = []
slist = []
tmp = []
hlist = []

if csvfile != None:
    n = 0
    for row in csv.reader(csvfile):
        for col in row:
            Lt = col.sprit().split('\t')
            if n == 0:
                tmp = Lt
                hlist = tmp[1:]
            else:
                cols = Lt
                slist.append(cols[0])
                dlist.append([float(str) for str in cols[1:]])
            n += 1
csvfile.close()

HN = vstack(hlist)
SN = vstack(slist)
R = vstack(dlist)
print('--分析対象となるデータ行列:R--')
print(R)

R_num_row = len(R)
R_num_col = size(R[0])
S = [0 for i in range(R_num_row)]
C = [0 for i in range(R_num_col)]

print('--カテゴリごとの点数が対角にある行列の逆行列:CT--')
for m in range(0, R_num_col):
    sum = 0
    for n in range(0, R_num_row):
        sum += R[n,m]
    C[m] = sum
CI = mat(diag(C)).I
print(CI)

print('--サンプルごとの点数の平方根が対角にある行列の逆行列:SI--')
for i in range(0, R_num_row):
    sum = 0
    for m in range(0, R_num_col):
        sum += R[n,m]
    S[n] = sqrt(sum)

SI = mat(diag(5)).I
print(SI)

print('--データ行列Rの転置行列:RT--')
RT = R.T
print(RT)

print('--SI*R*CI*RT*SI*v=r2*vを満たす固有値・固有ベクトルを求める--')
X = SI*R*CI*RT*SI
v, r2, vh = linalg.svd(X)
print('固有値:r2')
print(r2)

print('寄与率')
r2 = r2[1:]
r3 = r2
for n in range(0, size(r3)):
    r3[n] = r3[n] / r3.sum()
print(r3)

print('固有ベクトル:v')
print(v)

print('単相関係数:r')
r = sqrt(r2)
print(r)

print('--サンプルスコア'--)
x1 = sqrt(R.sum())*v[:,1:]
x1 = SI*x1
print(x1)

print('カテゴリスコア')
y1 = CI*R.T*x1
y2 = zeros((len(y1), size(y1[0])))
for j in range(0, size(y1[0])):
    y[:,j:j+1] = y1[:,j] / r[j]
y2 = mat(y2)
print(y2)

print('--散布図表示 (第一成分X軸・第二成分Y軸) --')
s_score_x = zeros((1,1))
s_score_y = zeros((1,1))
c_score_x = zeros((1,1))
c_score_y = zeros((1,1))

s_score_x = array(x1[:,0])
s_score_y = array(x1[:,1])
c_score_x = array(y2[:,0])
c_score_y = array(y2[:,1])

fig1 = plt.figure()

ax = fig1.add_subplot(111, autoscale_on=False, xlim(-1,5), ylim=(-4,3))
l, = plt.plot([], [], 'r-')

fig1 = plt.scatter(s_score_x, s_score_y, marker='O', edgecolors='r', facecolor='None', s=40, alpha=1, label='Sample Score')
fig1 = plt.scatter(c_score_x, c_score_y, marker='D', edgecolors='b', facecolor='None', s=40, alpha=1, label='Category Score')

fp = FontProperties(fname=r'C:\WINDOWS\Fonts\mshothic.ttc')

for n in range(0, size(s_score_x)):
    SN[n]
    ax.annotate(
            slist[n],
            xy=(s_score_x[n]-0.05, s_score_y[n]),
            xytext=(s_score_x[n]-0.2, s_score_y[n]),
            arrowprops=dict(
                arrowstyle='-',
                facecolor='black',
                ),
            horizontalaligment='right', verticalaligment='center',
            fontproperties=fp
    )

for n in range(0, size(c_score_x)):
    ax.annotate(
            hlist[n],
            xy=(c_score_x[n]+0.05, c_score_y[n]),
            xytext=(c_score_x[n]+0.2, c_score_y[n]),
            arrowprops=dict(
                arrowstyle='-',
                facecolor='black',
            ),
        horizontalaligment='left', verticalaligment='center',
        fontproperties=fp
    )

plt.xlim(-3,3)
plt.ylim(-3,3)
plt.axhline()
plt.axvline()

plt.grid(True)
plt.legend(loc='best')
plt.xlabel('成分1', fontproperties=fp)
plt.ylabel('成分2', fontproperties=fp)
plt.title('数量化三類', fontproperties=fp)
plt.show()
