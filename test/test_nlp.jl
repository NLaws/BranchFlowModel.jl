@testset "NLP" begin
    dssfilepath = "data/ieee13/IEEE13Nodeckt_no_trfxs.dss"
    net = BranchFlowModel.CommonOPF.dss_to_Network(dssfilepath)

    net.Vbase = 2400
    net.Sbase = 1e3
    net.Zbase = net.Vbase^2/net.Sbase
    net.bounds.v_upper = 1.1
    net.bounds.v_lower = 0.9
    net.bounds.s_upper =  net.Sbase * 100
    net.bounds.s_lower = -net.Sbase

    m = Model(Ipopt.Optimizer)

    BranchFlowModel.add_bfm_variables(m, net)
    BranchFlowModel.constrain_bfm_nlp(m, net)

end