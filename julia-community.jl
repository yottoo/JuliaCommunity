"""
This program is to discover and plot the communities of a network by leiden algorithm
NOTE: the leiden algorithm is implemented by the python package leidenalg, so
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
Github: https://github.com/yottoo/JuliaCommunity
"""

module JuliaCommunity

using LightGraphs, SimpleWeightedGraphs, GraphPlot
using PyCall
using DataFrames, CSV
using Statistics, StatsBase, Random
using Parameters, ProgressMeter
using Gadfly, Cairo, Compose

export JuliaCommunityInstance, 
        discover_communities, 
        optimise_resolution, 
        save,
        plot_community,
        compute_cluster_coef

const leiden = pyimport("leidenalg")
const ig = pyimport("igraph")

@with_kw mutable struct PartitionMethod
    louvain::String = "louvain"
    CPM::String = "CPM"
    modularity::String = "modularity"
end


@with_kw mutable struct JuliaCommunityInstance
    task_series::String = ""

    methods::PartitionMethod = PartitionMethod()
    method::String = "CPM"    # CPM, Modularity (for leiden algorithm) and Louvain
    """
    NETWORK: the network data.
        from: id of the node from
        to: id of the node to
        weight: edge weight (if edge_weighted is false, this column could be ignored)
    NODES: the data for nodes, with ID required
        id: node id
    NODE_IMPORTANCES: the importances of vertices
        id: node id
        importance
    """    
    network::DataFrame = DataFrame()
    graph = nothing
    igraph = nothing
    nodes::DataFrame = DataFrame()
    max_node_id::Int = 1
    node_importance_field::String = "importance"
    node_label_field::String = "id"
    
    edge_weighted::Bool = true
    node_weighted::Bool = false
    is_directed::Bool = true
    γ::Float16 = 0.001    

    """
    COMMUNITIES: communities discovered.
        c: community id (start from 1)
        size: community size
        cluster_cof: clustering coefficent of the community
    MEMBERSHIPS: memberships indicating the nodes belonging to communities
        id: node id
        c: community id   
    """
    communities::DataFrame = DataFrame()
    membership_vector::Array{Int} = []
    memberships::DataFrame = DataFrame()
    n_community::Int = 0

    modularity::Float64 = 0
    quality::Float64 = 0
end


function JuliaCommunityInstance(network::DataFrame; 
        nodes::DataFrame=DataFrame(), 
        node_label_field::String="id",
        node_importance_field="importance", 
        edge_weighted::Bool=true, 
        node_weighted::Bool=false, 
        is_directed::Bool=true, 
        method::String="CPM", 
        to_summarise_graph::Bool=true,
        task_series::String="",)
    filter!(row -> row.from > 0 && row.to > 0, network)    
    jc = JuliaCommunityInstance()
    jc.task_series = replace("_" * replace(task_series, "-" => "_"), "__" => "_")
    set_method(jc, method=method)
    print("\nLoading the network and nodes data...")
    jc.network = network
    check_network(jc)
    jc.nodes = nodes
    check_nodes(jc)
    jc.max_node_id = maximum(jc.nodes.id)
    jc.node_importance_field = node_importance_field
    jc.node_label_field = node_label_field
    jc.edge_weighted = edge_weighted
    jc.node_weighted = node_weighted
    jc.is_directed = is_directed
    
    if edge_weighted
        filter!(row -> row.weight > 0, jc.network) 
    #= ============================================
    else if !is_directed
        jc.network = vcat(jc.network, jc.network)
        unique(jc.network)
    ============================================ =#
    end

    print("\n\tBuilding the graph from the network...")
    if jc.is_directed && jc.edge_weighted
        jc.graph = SimpleWeightedDiGraph(jc.network.from, jc.network.to, jc.network.weight)
        jc.igraph = ig.Graph(zip(jc.network.from .- 1, jc.network.to .- 1), directed=true, edge_attrs=Dict("weight" => jc.network.weight))
    elseif jc.is_directed        
        jc.graph = SimpleDiGraph(jc.max_node_id)
        for i in 1:nrow(jc.network) add_edge!(jc.graph, jc.network.from[i], jc.network.to[i]) end
        jc.igraph = ig.Graph(zip(jc.network.from .- 1, jc.network.to .- 1), directed=true)
    elseif jc.edge_weighted
        jc.graph = SimpleWeightedGraph(jc.network.from, jc.network.to, jc.network.weight)
        jc.igraph = ig.Graph(zip(jc.network.from .- 1, jc.network.to .- 1), directed=false, edge_attrs=Dict("weight" => jc.network.weight))
    else
        jc.graph = SimpleGraph(jc.max_node_id)
        for i in 1:nrow(jc.network) add_edge!(jc.graph, jc.network.from[i], jc.network.to[i]) end
        jc.igraph = ig.Graph(zip(jc.network.from .- 1, jc.network.to .- 1), directed=false)
    end
      
    if to_summarise_graph summarise_graph(jc) end

    jc        
end


function check_network(jc::JuliaCommunityInstance)
    columns = names(jc.network)
    if isnothing(findfirst(name -> name == "from", columns)) throw(error("The network data must have a 'from' column.")) end
    if isnothing(findfirst(name -> name == "to", columns)) throw(error("The network data must have a 'to' column.")) end
    if jc.edge_weighted && isnothing(findfirst(name -> name == "weight", columns)) throw(error("The network data must have a 'weight' column.")) end    
end


function check_nodes(jc::JuliaCommunityInstance)
    columns = names(jc.nodes)
    if isnothing(findfirst(name -> name == "id", columns)) throw(error("The nodes data must have a 'id' column.")) end    
end

function set_method(jc::JuliaCommunityInstance; method::String="CPM")
    jc.method = method
    if (jc.method != jc.methods.louvain && jc.method != jc.methods.CPM && jc.method != jc.methods.modularity) 
        throw(error("The partition method has to be louvain, CMP or modularity.")) 
    end
end


function summarise_graph(jc::JuliaCommunityInstance)
    print("\n\tSummary of the graph built on the network: ")
    print("\n\t\tNodes: $(nv(jc.graph))[$(nrow(jc.nodes))]")
    print("\n\t\tEdges: $(ne(jc.graph))")
    print("\n\t\tDensity: $(density(jc.graph))")
    print("\n\t\tDiameter: $(jc.igraph.diameter())")
    print("\n\t\tAverage path length: $(jc.igraph.average_path_length(directed=jc.is_directed))")
    print("\n\t\tClustering coefficient: $(global_clustering_coefficient(jc.graph))")    
    print("\n\t\tModularity: $(modularity(jc.graph, fill(1, nv(jc.graph)), γ=1.0))\n")
end


function plot_network(jc::JuliaCommunityInstance; fig_path::String="fig", mute::Bool=false, node_size_smoother::Float64=0.5, edge_width_smoother::Float64=0.5, line_type::String="straight", arrow_length_frac::Float64=0.015, arrow_angle::Float64=π / 8)
    if !mute print("\n\nPlotting the network graph......\n\t")    end
    
    g = jc.graph

    if nv(g) > 5000 throw(error("It makes no sense to plot a huge network with more than 5000 vertices.")) end

    labels = jc.nodes[:, jc.node_label_field]
    nodesize = jc.node_weighted ? jc.nodes[:, jc.node_importance_field].^node_size_smoother : fill(1, nv(g))
    # labels[findall(size -> size < ceil(ne(jc.graph) / 100), nodesize)] .= ""
    # nodesize = log.(nodesize .+ 1)
    # edge_weights = jc.network[:weight].^0.72
    edge_weight = jc.edge_weighted ? jc.network.weight.^edge_width_smoother : fill(1, ne(g))
    # g = Graph(adjacency_matrix(g))
    node_colorss = ["orange", "purple", "turquoise", "green", "red", "blue", "violet", "olive", "tan", "magenta", "cyan", "pink", "gold"]
    edge_colors = ["navajowhite3", "coral2", "orange3", "yellow3", "yellowgreen", "turquoise", "lightskyblue", "mediumpurple1", "hotpink2", "tan3", "grey64"]
    shuffle!(node_colorss)
    # node_colors = colors[1 .+ Int.((length(colors) - 1) .* ceil.((nodesize .-  minimum(nodesize)) ./ (maximum(nodesize) - minimum(nodesize))))]
    # edge_colors = colors[1 .+ Int.((length(colors) - 1) .* ceil.((weights .-  minimum(weights)) ./ (maximum(weights) - minimum(weights))))]
    
    #= ========================================================================
    Further partition the target community and render each sub partition with a same color
    ======================================================================== =#
    igraph = jc.igraph
    
    categories = leiden.find_partition(igraph, leiden.ModularityVertexPartition)
    # categories = leiden.find_partition(igraph, leiden.CPMVertexPartition, resolution_parameter= 1 / nv(g))
    node_colors = node_colorss[(categories.membership .+ 1) .% length(node_colorss) .+ 1]
    # layout = (args...) -> spring_layout(args...; C = 12)   # where C influences the desired distance between nodes.

    run_label = "$(jc.method)-$(jc.γ)$(jc.task_series)"

    plot = nothing

    if jc.is_directed
        plot = gplot(g, nodesize=nodesize, nodelabel=labels, 
                    nodelabeldist=0.2, nodelabelangleoffset=π / 4,
                    nodelabelsize=nodesize, edgelinewidth=edge_weight, 
                    nodefillc=node_colors,
                    nodelabelc=node_colors, 
                    linetype=line_type,
                    edgestrokec=edge_colors[rand(1:length(edge_colors), 1)[1]], arrowlengthfrac=arrow_length_frac, arrowangleoffset=arrow_angle
                    ); 
                    # or linetype = "curve" or "straight" 
                    # edgestrokec = edge_colors,  
                    # edgestrokec=colors[rand(1:length(colors), 1)[1]]
    else
        plot = gplot(g, nodesize=nodesize, nodelabel=labels, 
                    nodelabeldist=0.2, nodelabelangleoffset=π / 4,
                    nodelabelsize=nodesize, edgelinewidth=edge_weight,
                    nodefillc=node_colors, 
                    nodelabelc=node_colors, 
                    linetype=line_type,
                    edgestrokec=edge_colors[rand(1:length(edge_colors), 1)[1]]);
    end
    if !ispath(fig_path) mkpath(fig_path) end
    filename = "$fig_path/network-graph-$run_label.svg"
    draw(SVG(filename, 20cm, 16cm), plot);
    open_file(filename)
end


"""
caculate_centralities:  Caculate the centralities of the network.
    Note: for weighted networks, eigenvector centrality measure may not lead to convergency.
"""
function caculate_centralities(jc::JuliaCommunityInstance; to_save::Bool=true) 
    print("\n\tCaculating the centralities of the target network......\n") 

    centralities = DataFrame(id=collect(vertices(jc.graph)))
    insertcols!(centralities, "btw_cent"  => betweenness_centrality(jc.graph),  
        "close_cent"  => closeness_centrality(jc.graph), 
        "in_cent"  => indegree_centrality(jc.graph), 
        "out_cent"  => outdegree_centrality(jc.graph),
        "eigen_cent"  => eigenvector_centrality(jc.graph),
        "katz_cent"  => katz_centrality(jc.graph, 0.3),
        "stress_cent"  => stress_centrality(jc.graph),
        "radia_cent"  => radiality_centrality(jc.graph))

    insertcols!(centralities, "btw_cent_rank" => competerank(centralities[:, "btw_cent"], rev=true),  
        "close_cent_rank" => competerank(centralities[:, "close_cent"], rev=true),
        "in_cent_rank" => competerank(centralities[:, "in_cent"], rev=true),
        "out_cent_rank" => competerank(centralities[:, "out_cent"], rev=true),
        "eigen_cent_rank" => competerank(centralities[:, "eigen_cent"], rev=true),
        "katz_cent"  * "_rank" => competerank(centralities[:, "katz_cent"], rev=true),
        "stress_cent_rank" => competerank(centralities[:, "stress_cent"], rev=true),
        "radia_cent"  * "_rank" => competerank(centralities[:, "radia_cent"], rev=true))

    centralities = leftjoin(jc.nodes, centralities, on=:id)
    
    if to_save
        save_path = "data/centralities$(jc.task_series).csv"
        centralities_merged =  ispath(save_path) ? hcat(CSV.File(save_path) |> DataFrame, centralities) : centralities
        CSV.write(save_path, centralities_merged)
    end

    centralities   
end


"""
_louvain:    the louvain algorithm implemented by leiden package    
"""
function _louvain(g)
    optimiser = leiden.Optimiser()
    partitions = leiden.ModularityVertexPartition(g)
    partitions_agg = partitions.aggregate_partition()
    while optimiser.move_nodes(partitions) > 0
        partitions.from_coarse_partition(partitions_agg)
        partitions_agg = partitions_agg.aggregate_partition()
    end
    partitions
end


"""
discover_communities:    discover the communities by leiden algorithm (both CMP and modularity method) 
    and louvain algorithm.
"""
function discover_communities(jc::JuliaCommunityInstance; mute::Bool=false)
    if !mute print("\nDiscovering the communities for the built network......") end
    
    partitions = nothing
    if jc.method == jc.methods.CPM
        partitions = leiden.find_partition(jc.igraph, leiden.CPMVertexPartition, resolution_parameter=jc.γ)
    elseif jc.method == jc.methods.modularity
        partitions = leiden.find_partition(jc.igraph, leiden.ModularityVertexPartition)
    elseif jc.method == jc.methods.louvain
        partitions = _louvain(jc.igraph)
    end
    
    if isnothing(partitions) return end

    jc.n_community = length(partitions)
    if !mute println("\t\t$(jc.n_community) communities have been discovered.") end
    
    # partitions = leiden.find_partition(g, leiden.ModularityVertexPartition, resolution_parameter=0.2)
    # The following code could not be done to Leiden community members 
    # for i in 1:length(partitions) partitions[i] .+= 1 end
    # partitions.membership .+= 1
    jc.membership_vector = partitions.membership .+ 1
    # println(partitions.membership)
    
    jc.modularity = modularity(jc.graph, jc.membership_vector, γ=1.0)
    jc.quality =  jc.method == jc.methods.CPM ? partitions.quality() : jc.modularity
    
    jc.memberships = DataFrame(id=1:maximum(vcat(jc.network.from, jc.network.to)), c=jc.membership_vector)
    # jc.memberships = DataFrame(id=1:nv(g), c=jc.membership_vector)
    jc.communities = DataFrame(c=1:jc.n_community, size=length.(partitions))
end

function optimise_resolution(jc::JuliaCommunityInstance; γ_from::Float64=0.0001, γ_end::Float64=0.01, γ_step::Float64=0.0001)
    println("\n")
    qualities = DataFrame(resolution=[], n_community=[], modularity=[], quality=[])
    jc_copy = deepcopy(jc)
    progress = Progress(length(γ_from:γ_step:γ_end), desc="Finding the best resolution γ for the Leiden-based community discovery algorithm: ")
    for γ  in γ_from:γ_step:γ_end
        jc_copy.γ = γ
        discover_communities(jc_copy, mute = true)
        push!(qualities, (γ, jc_copy.n_community, jc_copy.modularity, jc_copy.quality))
        next!(progress)
        #print("\n\t\tResolution: $γ: $(jc_copy.n_community) commxunities discovered; Modularity: $(jc_copy.modularity); CPM Quality: $(jc_copy.quality).") 
    end
    CSV.write("data/community_discover_optimisation-$(jc.method)$(jc.task_series).csv", qualities)

    p_modularity = plot(
        layer(qualities, x=:resolution, y=:modularity, Geom.line, Geom.point, Theme(default_color="blue")),        
        Guide.xticks(ticks=γ_from:γ_step * 2:γ_end),
        Guide.xlabel("resolution γ"),
        Guide.ylabel("modularity"),
        Theme(major_label_font_size=10pt),
        Scale.y_continuous(format=:plain)
    )
    fig_modularity = "fig/community_discover_optimisation-modularities-$(jc.method)$(jc.task_series).svg"
    draw(SVG(fig_modularity, 24cm, 16cm), p_modularity);
    open_file(fig_modularity)

    p_quality = plot(
        layer(qualities, x=:resolution, y=:quality, Geom.line, Geom.point, Theme(default_color="blue")),        
        Guide.xticks(ticks=γ_from:γ_step * 2:γ_end),
        Guide.xlabel("resolution γ"),
        Guide.ylabel("quality"),
        Theme(major_label_font_size=10pt),
        Scale.y_continuous(format=:plain)
    )
    fig_quality = "fig/community_discover_optimisation-qualities-$(jc.method)$(jc.task_series).svg"
    draw(SVG(fig_quality, 24cm, 16cm), p_quality);
    open_file(fig_quality)

    p_n_communities = plot(
        layer(qualities, x=:resolution, y=:n_community, Geom.line, Geom.point, Theme(default_color="blue")),        
        Guide.xticks(ticks=γ_from:γ_step * 2:γ_end),
        Guide.xlabel("resolution γ"),
        Guide.ylabel("number of communities"),
        Theme(major_label_font_size=10pt),
        Scale.y_continuous(format=:plain)
    )
    fig_n_communities = "fig/community_discover_optimisation-n-communities-$(jc.method)$(jc.task_series).svg"
    draw(SVG(fig_n_communities, 24cm, 16cm), p_n_communities);
    open_file(fig_n_communities)
end


function save(jc::JuliaCommunityInstance; file::String="data/communities")    
    print("\n\nSaving the communities discovery outcomes......\n")
    run_label = "$(jc.method)-$(jc.γ)$(jc.task_series)"
    CSV.write("$file-$run_label.csv", jc.communities)
    CSV.write("$file-memberships-$run_label.csv", jc.memberships)
end


function build_community_graph(jc::JuliaCommunityInstance, c::Int; to_save_data::Bool=true, mute::Bool=false)
    if !mute print("\tBuilding the graph for community $c......\n\t") end
    community = select(filter(row -> row.c == c, jc.memberships), Not(:c))
    if nrow(community) < 2 return (graph = nothing, edges = 1) end
    sort!(community, :id)
    community_size = nrow(community)

    nodes = DataFrame(_label=community.id, new_id=1:nrow(community))
    if jc.node_weighted
        nodes = select(innerjoin(nodes, select(jc.nodes, "id", jc.node_importance_field), on=:_label => :id), "_label", "new_id", jc.node_importance_field => "importance")
    end
    
    if jc.edge_weighted
        sub_network = rename(select(innerjoin(nodes, jc.network, on=:_label => :from), [:new_id, :to, :weight]), :new_id => "from")
        sub_network = rename(select(innerjoin(nodes, sub_network, on=:_label => :to), [:new_id, :from, :weight]), :new_id => "to")
    else
        sub_network = rename(select(innerjoin(nodes, jc.network, on=:_label => :from), [:new_id, :to]), :new_id => "from")
        sub_network = rename(select(innerjoin(nodes, sub_network, on=:_label => :to), [:new_id, :from]), :new_id => "to")
    end

    # There may be communities in which all nodes are not connected to others
    if nrow(sub_network) <= 0 return (graph = nothing, edges = 0) end

    if jc.node_label_field != "id"
        nodes = innerjoin(select(jc.nodes, :id, jc.node_label_field), nodes, on=:id => :_label)
        rename!(nodes, jc.node_label_field => "label")
    else
        rename!(nodes, :_label => "label")
    end

    graph = nothing
    if jc.is_directed && jc.edge_weighted
        graph = SimpleWeightedDiGraph(sub_network.from, sub_network.to, sub_network.weight)
    elseif jc.is_directed
        graph = SimpleDiGraph(community_size)
        for i in 1:nrow(sub_network) add_edge!(graph, sub_network.from[i], sub_network.to[i]) end
    elseif jc.edge_weighted
        graph = SimpleWeightedGraph(sub_network.from, sub_network.to, sub_network.weight)
    else
        graph = SimpleGraph(community_size)
        for i in 1:nrow(sub_network) add_edge!(graph, sub_network.from[i], sub_network.to[i]) end
    end

    if to_save_data
        CSV.write("data/graph_network_$(jc.method)-$(jc.γ)$(jc.task_series)-$c.csv", sub_network)
        if jc.node_weighted CSV.write("data/graph_nodes_$(jc.method)-$(jc.γ)$(jc.task_series)-$c.csv", nodes) end
    end

    if jc.node_weighted
        (graph = graph, node_labels = nodes.label, network = sub_network, community_size = community_size, node_weights = nodes.importance)
    else
        (graph = graph, node_labels = nodes.label, network = sub_network, community_size = community_size)
    end   
end

function plot_community(jc::JuliaCommunityInstance, c::Int; fig_path::String="fig", to_save_data::Bool=true, mute::Bool=false, node_size_smoother::Float64=0.5, edge_width_smoother::Float64=0.5, line_type="straight", arrow_length_frac::Float64=0.012, arrow_angle::Float64=π / 18)
    if !mute print("\n\nPlotting the community $c......\n\t")    end
    
    community_graph = build_community_graph(jc, c, to_save_data=to_save_data, mute=mute)
    if isnothing(community_graph) return end
    g = community_graph.graph
    if isnothing(g) return  end
    labels = community_graph.node_labels
    nodesize = jc.node_weighted ? community_graph.node_weights.^node_size_smoother : fill(1, nv(g))
    # labels[findall(size -> size < ceil(community.edges / 100), nodesize)] .= ""
    # nodesize = log.(nodesize .+ 1)
    # edge_weights = log.(community.network[:co_f] .+ 1)
    # edge_weights = community_graph.network[:weight].^0.72
    edge_weight = jc.edge_weighted ? community_graph.network.weight.^edge_width_smoother : fill(1, ne(g))
    # g = Graph(adjacency_matrix(g))
    node_colorss = ["orange", "purple", "turquoise", "green", "red", "blue", "violet", "olive", "tan", "magenta", "cyan", "pink", "gold"]
    edge_colors = ["navajowhite3", "coral2", "orange3", "yellow3", "yellowgreen", "turquoise", "lightskyblue", "mediumpurple1", "hotpink2", "tan3", "grey64"]
    shuffle!(node_colorss)
    # node_colors = colors[1 .+ Int.((length(colors) - 1) .* ceil.((nodesize .-  minimum(nodesize)) ./ (maximum(nodesize) - minimum(nodesize))))]
    # edge_colors = colors[1 .+ Int.((length(colors) - 1) .* ceil.((weights .-  minimum(weights)) ./ (maximum(weights) - minimum(weights))))]
    
    #= ========================================================================
    Further partition the target community and render each sub partition with a same color
    ======================================================================== =#
    network = community_graph.network
    igraph = nothing
    if jc.is_directed && jc.edge_weighted
        igraph = ig.Graph(zip(network.from .- 1, network.to .- 1), directed=true, edge_attrs=Dict("weight" => network.weight))
    elseif jc.is_directed        
        igraph = ig.Graph(zip(network.from .- 1, network.to .- 1), directed=true)
    elseif jc.edge_weighted
        igraph = ig.Graph(zip(network.from .- 1, network.to .- 1), directed=false, edge_attrs=Dict("weight" => network.weight))
    else
        igraph = ig.Graph(zip(network.from .- 1, network.to .- 1), directed=false)
    end

    categories = leiden.find_partition(igraph, leiden.ModularityVertexPartition)
    # categories = leiden.find_partition(igraph, leiden.CPMVertexPartition, resolution_parameter= 1 / nv(g))
    node_colors = node_colorss[(categories.membership .+ 1) .% length(node_colorss) .+ 1]  
    # layout = (args...) -> spring_layout(args...; C = 12)   # where C influences the desired distance between nodes.

    run_label = "$(jc.method)-$(jc.γ)$(jc.task_series)"

    plot = nothing

    if jc.is_directed
        plot = gplot(g, nodesize=nodesize, nodelabel=labels, 
                    nodelabeldist=0.2, nodelabelangleoffset=π / 4,
                    nodelabelsize=nodesize, edgelinewidth=edge_weight, 
                    nodefillc=node_colors, 
                    nodelabelc=node_colors, 
                    linetype=line_type,
                    edgestrokec=edge_colors[rand(1:length(edge_colors))], arrowlengthfrac=arrow_length_frac, arrowangleoffset=arrow_angle
                    ); 
                    # or linetype = "curve" or "straight" 
                    # edgestrokec = edge_colors,  
                    # edgestrokec=colors[rand(1:length(colors), 1)[1]]
    else
        plot = gplot(g, nodesize=nodesize, nodelabel=labels, 
                    nodelabeldist=0.2, nodelabelangleoffset=π / 4,
                    nodelabelsize=nodesize, edgelinewidth=edge_weight, 
                    nodefillc=node_colors, 
                    nodelabelc=node_colors, 
                    linetype=line_type,
                    edgestrokec=edge_colors[rand(1:length(edge_colors))]);
    end
    save_path = "$fig_path/community-$run_label" 
    if !ispath(save_path) mkpath(save_path) end
    draw(SVG("$save_path/$c.svg", 20cm, 16cm), plot);    
end


function compute_cluster_coef(jc::JuliaCommunityInstance)
    print("\nComputing the clustering coefficients for the partitioned communities......\n\t") 
    
    # If multiple threads are used to speed the caculation, then we have to bind the thread id with the corresponding index of the array.    
    cluster_cofs = zeros(jc.n_community)
    
    Threads.@threads for c in 1:jc.n_community
        community = build_community_graph(jc, c, to_save_data=false, mute=true)
        g = community.graph
        if isnothing(g)
            cluster_cofs[c] = community.edges
        else           
            cluster_cofs[c] = global_clustering_coefficient(g)
        end
    end

    insertcols!(jc.communities, :cluster_cof => cluster_cofs)
    avg_cluster_cof = mean(cluster_cofs)
    
    print("\tThe total clustering coefficient of the network  all the communities is $(global_clustering_coefficient(jc.graph)); The average clustering coefficient for all the communities is $avg_cluster_cof.")
end

function open_file(filename)
    if Sys.isapple()
        run(`open $(filename)`)
    elseif Sys.islinux() || Sys.isbsd()
        run(`xdg-open $(filename)`)
    elseif Sys.iswindows()
        run(`$(ENV["COMSPEC"]) /c start $(filename)`)
    else
        @warn "Showing plots is not supported on OS $(string(Sys.KERNEL))"
    end
end

end