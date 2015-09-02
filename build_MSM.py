import utilities as util


parent_dir = '/scratch/PI/rondror/DesRes-Simulations/b2AR_PNAS_2011a_Dror_Desres/DESRES-Trajectory_pnas2011a-A-no-water-no-lipid.1'
A_run_dir = '/DESRES-Trajectory_pnas2011a-A-no-water-no-lipid/pnas2011a-A-1-no-water-no-lipid'
dir_top = '/home/shenglan/topologies'

trajs_A_1 = util.load_trajs(A_run_dir,parent_dir,dir_top,load_stride = 2000)

