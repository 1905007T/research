import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#table_data = pd.read_table('ex100.txt')
#table_data

import csv
import matplotlib.pyplot as plt

for i in range(30,5000,30):
  # CSVファイルを読み込む
  with open(f'./graph/{i}.csv', 'r') as file:
      reader = csv.reader(file)
      # 最初の行をスキップする
      next(reader)
      # xとyの値をリストに格納する
      x_values = []
      y_values = []
      for row in reader:
          x_values.append(float(row[0]))
          y_values.append(float(row[1]))
  
  # グラフを描画する
  plt.plot(x_values, y_values)
  plt.xlabel('Cell Index')
  plt.ylabel('Electric Field(V/m)')
  plt.ylim(-70000,70000)
  
  # pngファイルとして保存する
  plt.savefig(f'./graph/{i}.png')
  plt.close()
  x_values.clear()
  y_values.clear()
  
import os
import imageio
import matplotlib.pyplot as plt

# アニメーションにする画像ファイルのリストを取得する
folder_path = './graph/'
files = sorted(os.listdir(folder_path), key=lambda x: int(x.split('.')[0]))
files = [folder_path + f for f in files if f.endswith('.png')]

# 画像を読み込んで、gifファイルとして保存する
with imageio.get_writer('animation.gif', mode='I', duration=1e-100) as writer:
    for filename in files:
        image = imageio.imread(filename)
        writer.append_data(image)
