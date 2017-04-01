# coding: UTF-8
#-------------------------------------------------------------------------------
# Name:        数量化3類実行プログラム
# Purpose:     数量化3類をPython Numpyで実行する
# Reference:   Excelで学ぶコレスポンデンス分析,高橋信,オーム社,2005
# Author:      Takao Aoki
# Created:     2013/01/20
# Copyright:   Takao Aoki(c)  2013
# Licence:      GNU GPL
# Description: 数量化3類とは、生データの回答者および変数に生データの情報が
# 十分に引き出された値を点グラフ化して表現する分析手法のこと
# 回答内容をダミー変数（０か１か）に置き換えて分析する。
#-------------------------------------------------------------------------------
from numpy import *
import numpy as np
from matplotlib.font_manager import FontProperties
from matplotlib.pyplot import figure, show
import matplotlib.pyplot as plt
import csv

# タブ区切りデータ読み込み
csvfile = open('data.txt')#分析対象データファイル
dlist = []#データ2次元リスト
slist = []#サンプル名称リスト
tmp =  []
hlist =  []#ヘッダーリスト

if csvfile != None:
    n = 0
    for row in csv.reader(csvfile):
        for col in row:
            Lt = col.decode('UTF-8').strip().split('\t')
            if n == 0:
                tmp = Lt
                hlist = tmp[1:]
            else:
                cols = Lt
                slist.append(cols[0])
                dlist.append([float(str) for str in cols[1:]])
            n += 1
csvfile.close
HN = vstack(hlist)
SN = vstack(slist)
R = vstack(dlist)
print u"--分析対象となるデータ行列:R--"
print R

R_num_row = len(R)#サンプルの数＝行数
R_num_col = size(R[0])#カテゴリの数＝列数
S = [0 for i in range(R_num_row)]
C = [0 for i in range(R_num_col)]

print u"--カテゴリごとの点数が対角にある行列の逆行列:CI--"
for m in range(0,R_num_col):
  sum = 0
  for n in range(0,R_num_row):
    sum += R[n,m]
  C[m] = sum
CI = mat(diag(C)).I#対角化・逆行列
print CI

print u"--サンプルごとの点数の平方根が対角にある行列の逆行列:SI--"
for n in range(0,R_num_row):
  sum = 0
  for m in range(0,R_num_col):
    sum += R[n,m]
  S[n] = sqrt(sum)

SI = mat(diag(S)).I#対角化・逆行列
print SI

print u"--データ行列Rの転置行列:RT--"
RT = R.T
print RT

print u"--SI*R*CI*RT*SI*v=r2*vを満たす固有値・固有ベクトルを求める--"
X = SI*R*CI*RT*SI
v,r2,vh = linalg.svd(X) #特異値分解のメソッドを使用する
print u"固有値:r2"
print r2

print u"寄与率"
r2 = r2[1:]#固有値1をスライスして削除
r3 = r2
for n in range(0,size(r3)):
  r3[n] = r3[n]/r3.sum()
print r3

print u"固有ベクトル:v"
print v

print u"単相関関係数:r"
r = sqrt(r2)
print r

print u"--サンプルスコア--"
x1 = sqrt(R.sum())*v[:,1:]#固有ベクトルに長さを乗する
x1 = SI*x1
print x1

print u"--カテゴリスコア--"
y1 = CI*R.T*x1
y2 = zeros((len(y1),size(y1[0])))
for j in range(0,size(y1[0])):#
  y2[:,j:j+1] = y1[:,j]/r[j]
y2 = mat(y2)
print y2

print u"--散布図表示（第1成分X軸・第2成分Y軸）--"
s_score_x = zeros((1,1))
s_score_y = zeros((1,1))
c_score_x = zeros((1,1))
c_score_y = zeros((1,1))

s_score_x = array(x1[:,0])#u"--サンプルスコアX軸--"
s_score_y = array(x1[:,1])#u"--サンプルスコアY軸--"
c_score_x = array(y2[:,0])#u"--カテゴリスコアX軸--"
c_score_y = array(y2[:,1])#u"--カテゴリスコアY軸--"

#グラフ表示
fig1=plt.figure()

#注記を表示
ax = fig1.add_subplot(111, autoscale_on=False, xlim=(-1,5), ylim=(-4,3))
l, = plt.plot([], [], 'r-')

#散布図のプロット
fig1 = plt.scatter(s_score_x,s_score_y, marker='o', edgecolors = 'r', facecolor='None', s=40, alpha = 1, label = "Sample Score")
fig1 = plt.scatter(c_score_x,c_score_y, marker='D', edgecolors = 'b', facecolor='None', s=40, alpha = 1, label = "Category Score")

#フォントはMSゴシック
fp = FontProperties(fname=r'C:\WINDOWS\Fonts\msgothic.ttc')

#サンプルスコアのラベル付
for n in range(0,size(s_score_x)):
    SN[n]
    ax.annotate(
            slist[n],
            xy=(s_score_x[n]-0.05, s_score_y[n]),
            xytext = (s_score_x[n]-0.2, s_score_y[n]),
            arrowprops=dict(
                arrowstyle="-",
                facecolor='black',
                ),
            horizontalalignment='right', verticalalignment='center',
            fontproperties=fp
    )
#カテゴリスコアのラベル付
for n in range(0,size(c_score_x)):
    ax.annotate(
            hlist[n],
            xy=(c_score_x[n]+0.05, c_score_y[n]),
            xytext = (c_score_x[n]+0.2, c_score_y[n]),
            arrowprops=dict(
                arrowstyle="-",
                facecolor='black',
                ),
            horizontalalignment='left', verticalalignment='center',
            fontproperties=fp
    )

#XY軸の最大値を決める処理

#座標軸を表示
plt.xlim(-3, 3)
plt.ylim(-3, 3)
plt.axhline()
plt.axvline()

#グラフ各種表示設定
plt.grid(True)
plt.legend(loc='best')
plt.xlabel(u'成分1', fontproperties=fp)
plt.ylabel(u'成分2', fontproperties=fp)
plt.title(u'数量化Ⅲ類', fontproperties=fp)
plt.show()


