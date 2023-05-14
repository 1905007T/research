import numpy as np
import matplotlib.pyplot as plt

# nodampenergy.csvファイルを読み込み、3つの列を配列として格納する
x, y, nodamp = np.loadtxt('nodampenergy.csv', delimiter=',', unpack=True)

# energy.csvファイルを読み込み、3つの列を配列として格納する
x, y, damp = np.loadtxt('energy.csv', delimiter=',', unpack=True)

decay=np.empty(len(nodamp))
print(type(nodamp))
print(type(decay))
for i in range(len(nodamp)):
  decay[i]=damp[i]/nodamp[i]
  print(i)
  print()
decay=10*np.log10(decay)
x=x/(2*np.pi)
y=y/(2*np.pi)

# カラーマップを作成する
plt.figure(figsize=(8, 6))
plt.scatter(x, y, c=decay, cmap='jet')
plt.colorbar()

# グラフタイトルや軸ラベルを設定する
plt.title('Reflection Ratio[dB]')
plt.xlabel('wpe')
plt.ylabel('wce')

# colormap.pngという名前で保存する
plt.savefig('colormap.png')
