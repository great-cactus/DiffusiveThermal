import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import Normalize
import glob
import time
import os; os.environ["OMP_NUM_THREADS"]='1'
plt.switch_backend('agg')
def main():
    target_file = 'data'
    files = glob.glob(target_file + '/*.vtk')
    files.sort()

    start = time.time()
    for file_name in files:
        print(file_name)
        tstp = int( file_name.split('/')[-1].split('.')[0][-6:])
        dt = 0.000002
        t = dt*tstp
        out_name = target_file + '/pngs/{:.4f}.png'.format(t)
        with open(file_name, 'r') as fr:
            line = fr.readlines()

        dimension = line[4].split()[1:]
        for i in range(len(dimension)): dimension[i] = int(dimension[i])

        points = dimension[0]*dimension[1]*dimension[2]
        grid = np.zeros( (points, 3) )

        for i, l in enumerate(line[6:6+points]):
            grid[i][0] = float(l.split()[0])
            grid[i][1] = float(l.split()[1])
            grid[i][2] = float(l.split()[2])

        data = np.zeros(points)
        for i, l in enumerate(line[6+points+3:]):
            data[i] = float(l)
        data = data.reshape( [dimension[0], dimension[1]] )

        x = np.linspace(-16, 16, dimension[0])
        y = np.linspace(-16, 16, dimension[1])

        fig, ax = plt.subplots( figsize=(10,10) )
        ax.pcolormesh(x, y, data, cmap='inferno', norm=Normalize(vmin=0, vmax=24))
        ax.set_title('t = {:.4f}'.format(t))
        plt.savefig(out_name, dpi=600)
        plt.close()

    print('visualization was finished. {:.4e}[sec] has passed.'.format(time.time()-start))

if __name__ == '__main__':
    main()
