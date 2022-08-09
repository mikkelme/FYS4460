from pylab import *
from scipy.sparse import spdiags, dia_matrix, coo_matrix
from scipy.sparse.linalg import spsolve
from scipy.ndimage import measurements
# from plot_set import *
# Written by Marin Soreng 2004
# Calculates the effective flow conductance Ceff of the
# lattice A as well as the potential V in every site .
def FIND_COND (A , X , Y ):
    V_in = 1.
    V_out = 0.
    # Calls MK_EQSYSTEM .
    B,C = MK_EQSYSTEM (A , X , Y )
    # Kirchhoff ’ s equations solve for V
    V = spsolve(B, C)
    # The pressure at the external sites is added
    # ( Boundary conditions )
    V = concatenate((V_in * ones (X), V, V_out * ones (X)))
    # Calculate Ceff
    # second-last X elements of V multiplied with second-last elem. of A
    # these are the second last column of the system
    # gives the conductivity of the system per row?
    Ceff = dot((V[-1-2*X:-1-X] - V_out).T, A[-1-2*X:-1-X, 1]) \
    / ( V_in - V_out )
    return V , Ceff

# Sets up Kirchoff ’ s equations for the 2 D lattice A .
# A has X * Y rows and 2 columns . The rows indicate the site ,
# the first column the bond perpendicular to the flow direction
# and the second column the bond parallel to the flow direction .
#
# The return values are [B , C ] where B * x = C . This is solved
# for the site pressure by x = B \ C .

def MK_EQSYSTEM (A , X , Y ):
    # Total no of internal lattice sites
    sites = X *( Y - 2)
    # Allocate space for the nonzero upper diagonals
    main_diag = zeros(sites)
    upper_diag1 = zeros(sites - 1)
    upper_diag2 = zeros(sites - X)
    # Calculates the nonzero upper diagonals
    main_diag = A[X:X*(Y-1), 0] + A[X:X*(Y-1), 1] + \
        A[0:X*(Y-2), 1] + A[X-1:X*(Y-1)-1, 0]

    upper_diag1 = A [X:X*(Y-1)-1, 0]
    upper_diag2 = A [X:X*(Y-2), 1]
    main_diag[where(main_diag == 0)] = 1
    # Constructing B which is symmetric , lower = upper diagonals .
    B = dia_matrix ((sites , sites)) # B *u = t
    B = - spdiags ( upper_diag1 , -1 , sites , sites )
    B = B + - spdiags ( upper_diag2 ,-X , sites , sites )
    B = B + B.T + spdiags ( main_diag , 0 , sites , sites )

    # Constructing C
    C = zeros(sites)
    # C = dia_matrix ( (sites , 1) )
    C[0:X] = A[0:X, 1]
    C[-1-X+1:-1] = 0*A [-1 -2*X + 1:-1-X, 1]

    return B , C

def sitetobond ( z ):
    # Function to convert the site network z(L,L) into a (L*L,2) bond
    # network
    # g [i,0] gives bond perpendicular to direction of flow
    # g [i,1] gives bond parallel to direction of flow
    # z [ nx , ny ] -> g [ nx * ny , 2]
    nx = size (z ,1 - 1)
    ny = size (z ,2 - 1)
    N = nx * ny
    gg_r = zeros ((nx , ny)) # First , find these
    gg_d = zeros ((nx , ny )) # First , find these
    gg_r [:, 0:ny - 1] = z [:, 0:ny - 1] * z [:, 1:ny]
    gg_r[:,ny -1]=z[:,ny -1]
    gg_d [0:nx - 1, :] = z [0:nx - 1, :] * z [1:nx, :]
    gg_d [nx - 1, :] = 0

    # Then , concatenate gg onto g
    g = zeros ((nx *ny ,2))
    g [:, 0] = gg_d.reshape(-1,order='F').T
    g [:, 1] = gg_r.reshape(-1,order='F').T



    return g

def coltomat (z, x, y):
    # Convert z(x*y) into a matrix of z(x,y)
    # Transform this onto a nx x ny lattice
    g = zeros ((x , y))
    for iy in range(1,y):
        i = (iy - 1) * x + 1
        ii = i + x - 1
        g[: , iy - 1] = z[ i - 1 : ii]
    return g


if __name__ == "__main__":
    # Generate spanning cluster (l - r spanning )
    lx = 400
    ly = 400
    p = 0.5927
    ncount = 0
    perc = []

    while (len(perc)==0):
        ncount = ncount + 1
        if (ncount >10):
            break
        np.random.seed(174412678) #Fixed seed <----------
        # np.random.seed(174412674) #Fixed seed <----------

        z=rand(lx,ly)<p

        lw,num = measurements.label(z)
        perc_x = intersect1d(lw[0,:],lw[-1,:])
        perc = perc_x[where(perc_x > 0)]
        print("Percolation attempt", ncount)
    zz = asarray((lw == perc[0]))
    # zz now contains the spanning cluster
    zzz = zz.T # Transpose
    g = sitetobond ( zzz ) # Generate bond lattice
    V, c_eff = FIND_COND (g, lx, ly) # Find conductivity
    x = coltomat ( V , lx , ly ) # Transform to nx x ny lattice
    V = x * zzz
    g1 = g[:,0]
    g2 = g[: ,1]
    z1 = coltomat( g1 , lx , ly )
    z2 = coltomat( g2 , lx , ly )


    # Calculate current from top to down from the potential
    f2 = zeros ( (lx , ly ))
    for iy in range(ly -1):
        f2[: , iy ] = ( V [: , iy ] - V [: , iy +1]) * z2 [: , iy ]
    # Calculate current from left to right from the potential
    f1 = zeros ( (lx , ly ))
    for ix in range(lx-1):
        f1[ ix ,:] = ( V [ ix ,:] - V [ ix +1 ,:]) * z1 [ ix ,:]
    # Find the sum of (absolute) currents in and out of each site
    fn = zeros (( lx , ly ))
    fn = fn + abs ( f1 )
    fn = fn + abs ( f2 )
    # Add for each column (except leftmost) the up-down current, but offset
    fn [: ,1: ly ] = fn [: ,1: ly ] + abs ( f2 [: ,0: ly -1])
    # For the left-most one, add the inverse potential
    # multiplied with the spanning cluster bool information
    fn [: ,0] = fn [: ,0] + abs (( V [: ,0] - 1.0)*( zzz [: ,0]))
    # For each row (except topmost) add the left-right current, but offset
    fn [1: lx ,:] = fn [1: lx ,:] + abs ( f1 [0: lx -1 ,:])



    # Singly connected
    zsc = fn > (fn.max() - 1e-6)
    # Backbone
    zbb = fn>1e-6
    # Combine visualizations
    ztt = ( zzz*1.0 + zsc*2.0 + zbb*3.0 )
    zbb = zbb / zbb.max()

    # cmap = plt.cm.OrRd

    cmap = plt.cm.viridis
    # cmap.set_under(color='#FFFFFF')
    # vmin = 1e-8
    # vmin = 0



    PowerPoint_plot = True
    Save = True
    if PowerPoint_plot:
        figsize = (6,6)

        figure(num = 0, figsize=figsize) #Spanning cluster
        imshow(zzz, interpolation='nearest', cmap=cmap)
        title("Spanning cluster")
        plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)
        if Save:
            plt.savefig("../article/figures/flow_solve_spanning_cluster.pdf", bbox_inches="tight")


        figure(num = 1, figsize=figsize) #Potential
        imshow(V, interpolation='nearest', cmap=cmap)
        title("Potential")
        plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)
        if Save:
            plt.savefig("../article/figures/flow_solve_potential.pdf", bbox_inches="tight")


        figure(num = 2, figsize=figsize) #Current
        imshow(fn, interpolation='nearest', cmap=cmap)
        title ("Current")
        plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)
        if Save:
            plt.savefig("../article/figures/flow_solve_current.pdf", bbox_inches="tight")


        figure(num = 3, figsize=figsize) #Current (lowered zero)
        fn_lowered = np.where(fn < 1e-8, fn-1/10*fn.max(), fn)

        imshow(fn_lowered, interpolation='nearest', cmap=cmap)
        title ("Current")
        plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)
        if Save:
            plt.savefig("../article/figures/flow_solve_current_lowered_zero.pdf", bbox_inches="tight")


        # figure(num = 3, figsize=figsize) #Current white
        # cmap.set_under(color='#FFFFFF')
        # vmin = 0.0
        # imshow(fn, interpolation='nearest', cmap=cmap, vmin = vmin)
        # title ("Current (BB)")
        # plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)
        # if Save:
        #     plt.savefig("../article/figures/flow_solve_current_white.pdf", bbox_inches="tight")


        figure(num = 4, figsize=figsize) #SC, BB, DE
        imshow(ztt, interpolation='nearest', cmap=cmap)
        title ("SC, BB and DE")
        plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)
        if Save:
            plt.savefig("../article/figures/flow_solve_SC_BB_DE.pdf", bbox_inches="tight")


        figure(num = 5, figsize=figsize) #BB, SC
        rgba = cmap(0.0)
        cmap.set_under(color=rgba)
        imshow(ztt, interpolation='nearest', cmap=cmap, vmin = 2)
        title ("BB, SC")
        plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)
        if Save:
            plt.savefig("../article/figures/flow_solve_BB_CC.pdf", bbox_inches="tight")




        figure(num = 6, figsize=figsize) #SC
        imshow(zsc, interpolation='nearest', cmap=cmap)
        title ("SC")
        plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)
        if Save:
            plt.savefig("../article/figures/flow_solve_SC.pdf", bbox_inches="tight")

        show()



    else:
        # Plot results
        figure(figsize=(6,6))

        ax = subplot(221)
        imshow(zzz, interpolation='nearest', cmap=cmap, vmin = vmin)
        title("Spanning cluster")

        subplot(222, sharex=ax, sharey=ax)
        imshow(V, interpolation='nearest', cmap=cmap, vmin = vmin)
        title("Potential")

        subplot(223, sharex=ax, sharey=ax)
        imshow(fn, interpolation='nearest', cmap=cmap, vmin = vmin)
        title ("Current (BB)")

        subplot(224, sharex=ax, sharey=ax)
        imshow(ztt, interpolation='nearest', cmap=cmap, vmin = vmin)
        title (" SC, BB and DE ")

        plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)
        # plt.savefig("../article/figures/flow_solve.pdf", bbox_inches="tight")
        show()
