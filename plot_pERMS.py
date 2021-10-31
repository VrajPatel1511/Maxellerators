import numpy as np
from matplotlib import animation
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.image as mgimg
import matplotlib.ticker
import matplotlib as mpl

fig, axs = plt.subplots()
plt.rcParams.update({'font.size': 16})
n = 0
m = 0
with open('./panimE/0.dat') as file:
    for line in file:
        m += 1 
print (m)
y = np.linspace(1, m, m)

with open('./panimE/0.dat') as file:
    line=file.readline()
    n=len(line.split())
#        n += 1 
print (n)
x = np.linspace(1, n, n)
X, Y = np.meshgrid(x, y)
#factor2=4*0.525e-5
factor2=1
def init():
    plt.xlabel('dist in lamb(m)')
    plt.ylabel('dist in lamb(m)')
    fmt = matplotlib.ticker.ScalarFormatter(useMathText=True)
    fmt.set_powerlimits((0,0))
    cbar = fig.colorbar(im,format=fmt)
    cbar.update_ticks()
    
ims = []
#ims2 = [1,2,3,4,5,6,7,8,9,10]
ims2 = [12]
#for i in range(27,28):
for i in ims2:
    with open('./panimE/' + str(i) + '.dat') as file:
        z = [[float(digit) for digit in line.split()] for line in file]

    maxz = np.amax(z)
#    maxz = 5.56E6
    minz = np.amin(z)
    levels = np.linspace(minz, maxz, 250)
    im = axs.contourf(factor2*X,factor2*Y, z, levels=levels,cmap="hot",alpha=0.85)
#    im = axs.contourf(factor2*X,factor2*Y, z, levels=levels)
    if i==12: 
   #    init()
   # plt.xlabel('dist in lamb(m)')
   # plt.ylabel('dist in lamb(m)')
   # fmt=matplotlib.ticker.ScalarFormatter(useMathText=True)
   # fmt.set_powerlimits((0,0))
   # cbar=fig.colorbar(im,format=fmt)
   # cbar.update_ticks()
     init()
    plt.savefig('./panimE/' + str(i) + '.png',dpi=750)
    
for p in range(1, 1):
    ## Read in picture
    fname = "./panimE/" + str(p)+ ".png"
    img = mgimg.imread(fname)
    imgplot = plt.imshow(img)

    # append AxesImage object to the list
    ims.append([imgplot])

## create an instance of animation
#ani = animation.ArtistAnimation(fig, ims, interval=400, blit=True, repeat_delay=1000)
#ani.save('/home/manyu/BTP/cuda_code/c_animE/animated_efield_p750.mp4')

plt.title('Electric Field with Plasma density 1e15')
#plt.xlabel('x grid dimension')
#plt.ylabel('y grid dimension')

#plt.show()
