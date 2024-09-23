import numpy as np
import math
import time
import matplotlib.pyplot as plt

'''
测试用数据集:
杭州市区经纬度:(120.20000,30.26667)     (作为起点)
上海市区经纬度:(121.43333,31.80000)   
南京市区经纬度:(118.78333,32.05000)
郑州市区经纬度:(113.65000,34.76667)
太原市区经纬度:(112.53333,37.86667)
石家庄市经纬度:(114.48333,38.03333)
天津市区经纬度:(117.20000,39.13333)
北京市区经纬度:(116.41667,39.91667)
沈阳市区经纬度:(123.38333,41.80000)
长春市区经纬度:(125.35000,43.88333)
哈尔滨市经纬度:(126.63333,45.75000)
'''

###Part 1   构建数据集
cn=int(input("请输入途径点个数（不含起点）："))          #city number 途径的城市数（不含起点）

##1_1 定义haversine函数_计算两点之间球面距离（经纬度作为参数）
def haversine(spot1,spot2):  # 经度1，纬度1，经度2，纬度2 （十进制度数）
    # 将十进制度数转化为弧度
    lon1, lat1, lon2, lat2 = map(math.radians, [spot1[0], spot1[1], spot2[0], spot2[1]])
    # haversine公式
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = math.sin(dlat / 2) ** 2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2) ** 2
    c = 2 * math.asin(math.sqrt(a))
    r = 6371  # 地球平均半径，单位为km
    distance=round(c * r * 1000,3)
    return distance

##1_2 手动输入起点+cn个途径点的经纬度，构建spot字典存储数据
spot={}
for i in range(cn+1):
    spot[i]=list(map(float,input("请输入第{}个点的经纬度（用半角逗号分割）:".format(i)).split(",")))
#显示输出
print("起点（第0点）和途经点的经纬度：")
print('\r')
for i in range (cn+1):
    print("第{}点：{},{}".format(i,spot[i][0],spot[i][1]))

##1_3 构造cn x cn矩阵，存储集合中两两点间距离
#构造矩阵并赋初值0
dis=np.empty((cn+1,cn+1))
dis.fill(0)
#赋实际距离
for i in range(cn+1):
    for j in range(i+1,cn+1):
        dis[i][j]=dis[j][i]=haversine(spot[i],spot[j])
        #dis[i][j] = dis[j][i] =int(input("请输入第{},{}点之间距离：".format(i,j)))  #10点共线测试用，此时需将1_1与1_2代码注释掉
#显示输出
print(dis)


start_time=time.time()     #计算运行时长
###Part 2   构建动态规划dp矩阵
cln=2**cn                  #dp矩阵column number （含有m个元素的集合，其子集数=2^(m-1)）
dp=np.empty((cn+1,cln))
dp.fill(float("inf"))
#初始化赋值dp[i][0]=第i点到起点的距离
for i in range (cn+1):
    dp[i][0]=dis[i][0]
#求解dp[i][j]，外循环更新列，内循环更新行
for j in range(1,cln):
    for i in range(cn+1):
        #如果集合j中包含第i点，则跳过此次循环
        if i!=0 and (j>>(i-1))&1==1 :
            continue
        for k in range(1,cn+1):
            if (j>>(k-1)&1)==0:
                continue
            elif dp[i][j] > dis[i][k]+dp[k][j^(1<<(k-1))] :
                dp[i][j] = dis[i][k] + dp[k][j ^ (1 << (k - 1))]
#显示输出
print(dp)


###Part 3  由dp表反推最短路径
#构建optimised path列表，依次保存经过点的序号
op_path=[]
#初始化变量
sec=cln-1           #初始集合包含所有途径点
pioneer=0           #上一个经过的点
node=0              #现在要去的点
min=float("inf")    #子问题路径的最小值
#循环求解各子问题的最优解
while sec!=0:
    #寻找现在要去的点，使子问题路径最短
    for i in range(1,cn+1):
        if (sec>>(i-1))&1 ==1:
            if min > dis[pioneer][i]+dp[i][sec^(1<<(i-1))] :
                min=dis[pioneer][i]+dp[i][sec^(1<<(i-1))]
                node=i
    #将该点存入optimised path
    op_path.append(node)
    #更新相关变量，准备下一次循环
    pioneer=node
    sec=sec^(1<<(node-1))
    min=float("inf")

end_time=time.time()             #计算运行时长
total_time=end_time-start_time
#显示输出
print("最短路径（最优解）为：")
print("起点-->",end='')
for i in range(len(op_path)):
    print("第{}点".format(op_path[i]),end="-->")
print("起点")
print("最短路径长度={}km".format(dp[0][cln-1]))
print("运行时长={}s".format(total_time))


###Part4 可视化图形输出(二维直角坐标近似画图，x轴为经度，y轴为纬度)
##4_1 设置figure
plt.figure()
plt.title("TSP_DP_Optimised Path")
plt.xlabel("Longitude")
plt.ylabel("Latitude")
##4_2 画出起点和途经点（起点☆，途经点●）
for i in range(cn+1):
    if i==0:
        plt.plot(spot[i][0],spot[i][1],marker="*",c='r')
    else:
        plt.plot(spot[i][0], spot[i][1], marker="o",c='r')
##4_3 依次画出路径
op_path.insert(0,0)
for i in range(len(op_path)-1):
    plt.plot((spot[op_path[i]][0], spot[op_path[i+1]][0]), (spot[op_path[i]][1], spot[op_path[i+1]][1]), linewidth=1, linestyle='-', c='b')
plt.plot((spot[op_path[-1]][0],spot[0][0]),(spot[op_path[-1]][1],spot[0][1]),linewidth=1,linestyle='-',c='b')
#输出最终图像
plt.show()
