import numpy as np
from matplotlib import animation
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.animation as animation
import matplotlib.image as mgimg

fig, axs = plt.subplots()
plt.rcParams.update({'font.size': 12})
n = 0
m = 0
with open('diffux1.dat') as file:
    for line in file:
        m += 1 
print(m)
y = np.linspace(1, m, m)

with open('diffux1.dat') as file:
    line=file.readline()
    n=len(line.split())
#        m += 1 
print(n)
x = np.linspace(1, n, n)
X, Y = np.meshgrid(x, y)
#factor2=4*0.525e-5
factor1=1/128
factor2=9.091e-12
def init():
    cbar = fig.colorbar(im)
#    plt.xlabel('dist in lamb(m)')
#    plt.ylabel('dist in lamb(m)')
#    norm=mpl.colors.Normalize(vmin=0.0E21,vmax=5.26E21)
#    cbar = fig.colorbar(im,ticks=[0.0E21,5.26E21])
#    cbar = fig.colorbar(im)
#    cbar.mappable.set_clim(vmin=0.0E21,vmax=5.26E21)
#    cbar.ax.set_yticklabels([0.0E21,4.00E21])
#    cbar.set_label('V/m')
    
ims = []
#ims2 = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
#ims2 = [1,2,3,4,5,6,7,8]
#for i in range(27,28):
#for i in ims2:
with open('diffux1.dat') as file:
     z = [[float(digit) for digit in line.split()] for line in file]

maxz = np.amax(z)
print(maxz)
#    maxz = 5.655E21
minz = np.amin(z)
print(minz)
levels = np.linspace(minz, maxz, 250)
im = axs.contourf(factor1*X,factor2*Y, z, levels=levels,cmap="jet",alpha=0.85)
#    im = axs.contourf(factor2*X,factor2*Y, z, levels=levels)
#    if i==8: 
init() 
    #init() 
   #     init()
plt.savefig('diffux1.png',dpi=750)
    
for p in range(1, 1):
    ## Read in picture
    fname = "diffux1.png"
    img = mgimg.imread(fname)
    imgplot = plt.imshow(img)

    # append AxesImage object to the list
    ims.append([imgplot])

## create an instance of animation
#ani = animation.ArtistAnimation(fig, ims, interval=400, blit=True, repeat_delay=1000)
#ani.save('/home/manyu/BTP/cuda_code/c_anim/animated_efield_p750.mp4')

plt.title('Electric Field with Plasma density 1e15')
#plt.xlabel('x grid dimension')
#plt.ylabel('y grid dimension')

#plt.show()
