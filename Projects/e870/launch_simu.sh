# # 12Be(p,a)9Li
# npsimulation -D DetectorConfiguration/MUGAST_LISE.detector -E reaction/12Bepalpha.reaction -O MugastAtLise12Bepalpha.root -B simu12Bepalpha.mac
# ninja && npanalysis -D DetectorConfiguration/MUGAST_LISE.detector -E reaction/12Bepalpha.reaction -T /Users/valerian/Software/nptool_gitlab/nptool/Outputs/Simulation/MugastAtLise12Bepalpha.root SimulatedTree -O AnaMugastAtLise12Bepalpha.root

# # 10Be(d,6Li)
# npsimulation -D DetectorConfiguration/MUGAST_LISE_CD2.detector -E reaction/10Bed6Li.reaction -O MugastAtLise10Bed6Li.root -B simu10Bed6Li.mac
# ninja && npanalysis -D DetectorConfiguration/MUGAST_LISE_CD2.detector -E reaction/10Bed6Li.reaction -T /Users/valerian/Software/nptool_gitlab/nptool/Outputs/Simulation/MugastAtLise10Bed6Li.root SimulatedTree -O AnaMugastAtLise10Bed6Li.root

# # 10Be(p,a)
# npsimulation -D DetectorConfiguration/MUGAST_LISE.detector -E reaction/10Bepalpha.reaction -O MugastAtLise10Bepalpha.root -B simu10Bepalpha.mac
# ninja && npanalysis -D DetectorConfiguration/MUGAST_LISE.detector -E reaction/10Bepalpha.reaction -T /Users/valerian/Software/nptool_gitlab/nptool/Outputs/Simulation/MugastAtLise10Bepalpha.root SimulatedTree -O AnaMugastAtLise10Bepalpha.root


# # 10Be(d,6Li)
# npsimulation -D DetectorConfiguration/MUGAST_LISE_CD2.detector -E reaction/10Bed6Li.reaction -O MugastAtLise10Bed6Li_thin_plus.root -B simu10Bed6Li_plus.mac
# ninja && npanalysis -D DetectorConfiguration/MUGAST_LISE_CD2.detector -E reaction/10Bed6Li.reaction -T /Users/valerian/Software/nptool_gitlab/nptool/Outputs/Simulation/MugastAtLise10Bed6Li_thin_plus.root SimulatedTree -O AnaMugastAtLise10Bed6Li_thin_plus.root

# 10Be(CD2)
# npsimulation -D DetectorConfiguration/MUGAST_LISE_CD2_thick.detector -E reaction/10Be.reaction -O MugastAtLise10Be_thick_plus.root -B simu10Bed6Li_plus.mac

# 10Be(d,6Li)
# npsimulation -D DetectorConfiguration/MUGAST_LISE_CD2_thick.detector -E reaction/10Be.reaction -O MugastAtLise10Be_thick_plus.root -B simu10Bed6Li_plus.mac
# npsimulation -D DetectorConfiguration/MUGAST_LISE_CD2_thick.detector -E reaction/10Bed6Li.reaction -O MugastAtLise10Bed6Li_thick_plus.root -B simu10Bed6Li_plus.mac

# 10Be(CH2)
# npsimulation -D DetectorConfiguration/MUGAST_LISE.detector -E reaction/10Be.reaction -O MugastAtLise10BeCH2.root -B simu10Bed6Li_plus.mac

# 10Be(p,alpha)
# npsimulation -D DetectorConfiguration/MUGAST_LISE.detector -E reaction/10Bepalpha.reaction -O MugastAtLise10Bepalpha.root -B simu10Bed6Li_plus.mac
# ninja && npanalysis -D DetectorConfiguration/MUGAST_LISE.detector -E reaction/10Bepalpha.reaction -T /Users/valerian/Software/nptool_gitlab/nptool/Outputs/Simulation/MugastAtLise10Bepalpha.root SimulatedTree -O AnaMugastAtLise10Bepalpha.root
# npanalysis -D DetectorConfiguration/MUGAST_LISE_CD2_thick.detector -E reaction/10Bed6Li.reaction -T /Users/valerian/Software/nptool_gitlab/nptool/Outputs/Simulation/MugastAtLise10Bed6Li_thick_plus.root SimulatedTree -O AnaMugastAtLise10Bed6Li_thick_plus.root
# npsimulation -D DetectorConfiguration/MUGAST_LISE.detector -E reaction/10Bepalpha.reaction -O MugastAtLise10Bepalpha.root -B simu10Bepalpha.mac

# npsimulation -D DetectorConfiguration/MUGAST_LISE_CD2_thick.detector -E reaction/10Bed6Li.reaction -O MugastAtLise10Bed6LiThick.root -B simu10Bed6Li.mac
# ninja && npanalysis -D DetectorConfiguration/MUGAST_LISE_CD2.detector -E reaction/10Bed6Li.reaction -T /Users/valerian/Software/nptool_gitlab/nptool/Outputs/Simulation/MugastAtLise10Bed6Li.root SimulatedTree -O AnaMugastAtLise10Bed6LiThick.root

ninja && npanalysis -D DetectorConfiguration/MUGAST_LISE_CD2.detector -E reaction/10Bed6Li.reaction -T /Users/valerian/Software/nptool_gitlab/nptool/Outputs/Simulation/MugastAtLise10Bed6Li.root SimulatedTree -O test.root

















