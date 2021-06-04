# JuliaCommunity
A julia wrapper of <a href='https://github.com/vtraag/leidenalg'>Leiden algorithm</a> to discover and plot the communities of a network.
<br>NOTE: the leiden algorithm is implemented by the python package leidenalg, so
        before doing community discovery, Conda, PyCall have to be installed as
        follows:

            import Pkg
            Pkg.add("Conda")
            Pkg.add("PyCall")
            Pkg.build("PyCall")
            
            using Conda

            Conda.pip_interop(true)
            #Conda.pip("install", "scipy")
            #Conda.pip("install", "numpy")
            Conda.pip("install", "leidenalg")

Contributors: Xiaoshan Nian
Date: August, 2020
Email: cen@njust.edu.cn
