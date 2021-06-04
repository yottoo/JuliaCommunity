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


<br><b>A Demo</b>
<pre>
include("julia-community/julia-community.jl")
using DataFrames

import .JuliaCommunity as juliac

nodes = DataFrame(id = [1,2,3,4,5,6,7,8,9,10], 
                label = ["a","b","c","d","e","f","g","h","i","j"], 
                importance=[1,5,6,5,5,4,3,3,2,2])
network = DataFrame(from = [1,2,2,2,3,3,4,4,4,5,5,5,6,6,8,8,10,10],
                      to = [3,3,4,9,2,4,2,3,7,3,6,7,5,8,5,6,7,9],
                  weight = [1,2,3,1,4,4,2,6,1,1,7,2,5,3,2,4,3,2])

jc = juliac.JuliaCommunityInstance(network, nodes = nodes, node_label_field = "label", node_weighted = true, to_summarise_graph = false, task_series = "test")
juliac.plot_network(jc, line_type="curve", node_size_smoother = 0.8, edge_width_smoother = 1.2)
jc.γ = 0.1
juliac.community_discover(jc)

println(jc.communities)
println(jc.memberships)

juliac.plot_community(jc, 1, line_type="curve")
</pre>
