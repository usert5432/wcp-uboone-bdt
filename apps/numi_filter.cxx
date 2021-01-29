// cz: code modified from tutorials/tmva/TMVAClassification.C

#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <set>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TMath.h"

#include "WCPLEEANA/tagger.h"

#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"


#include "WCPLEEANA/eval.h"

using namespace std;
using namespace LEEana;

#include "WCPLEEANA/bdt.h"
#include "WCPLEEANA/pot.h"
#include "WCPLEEANA/pfeval.h"
#include "WCPLEEANA/kine.h"
#include "WCPLEEANA/cuts.h"

int main( int argc, char** argv )
{
  if (argc < 4) {
    std::cout << "numi_filter #input_file #prefix_outfile -f[#filter_level] -r[#run_filter]" << std::endl;
    
    return -1;
  }

  TString input_file = argv[1];
  TString prefix_out = argv[2];

  Int_t filter_level = 1;
  Int_t run_filter = 0;
  
  for (Int_t i=1;i!=argc;i++){
    switch(argv[i][1]){
    case 'r':
      run_filter = atoi(&argv[i][2]);
      break;
    case 'f':
      filter_level = atoi(&argv[i][2]);
      break;
    }
  }
  TString outfile_name;


  std::vector<int> good_run_list_vec{4952, 4953, 4954, 4955, 4957, 4958, 4961, 4962, 4966, 4967, 4968, 4969, 4971, 4974, 4975, 4977, 4978, 4979, 4981, 4982, 4983, 4986, 4987,
4988, 4989, 4991, 4992, 4995, 4997, 4998, 4999, 5000, 5001, 5002, 5005, 5009, 5010, 5011, 5012, 5013, 5015, 5016, 5017, 5019, 5021, 5022, 5023, 5024,
5025, 5027, 5029, 5031, 5034, 5036, 5037, 5038, 5039, 5041, 5042, 5044, 5046, 5047, 5048, 5049, 5051, 5055, 5056, 5059, 5060, 5061, 5062, 5064, 5065,
5066, 5067, 5068, 5069, 5070, 5071, 5074, 5075, 5076, 5077, 5078, 5079, 5082, 5084, 5086, 5087, 5089, 5090, 5091, 5092, 5093, 5095, 5097, 5098, 5099,
5100, 5102, 5103, 5104, 5106, 5108, 5109, 5110, 5114, 5121, 5122, 5124, 5125, 5127, 5128, 5130, 5133, 5134, 5135, 5136, 5137, 5138, 5139, 5142, 5143,
5144, 5145, 5146, 5147, 5151, 5153, 5154, 5155, 5157, 5159, 5160, 5161, 5162, 5164, 5165, 5166, 5167, 5168, 5169, 5170, 5171, 5176, 5177, 5179, 5181,
5182, 5183, 5184, 5185, 5187, 5189, 5190, 5191, 5192, 5194, 5195, 5197, 5198, 5201, 5203, 5204, 5205, 5207, 5208, 5211, 5212, 5213, 5214, 5215, 5216,
5217, 5219, 5222, 5223, 5226, 5227, 5229, 5233, 5235, 5237, 5262, 5263, 5264, 5265, 5266, 5267, 5268, 5269, 5270, 5271, 5272, 5273, 5274, 5275, 5277,
5278, 5279, 5280, 5281, 5315, 5320, 5321, 5322, 5326, 5328, 5329, 5330, 5331, 5332, 5333, 5334, 5337, 5338, 5339, 5340, 5341, 5343, 5344, 5345, 5347,
5348, 5349, 5351, 5353, 5354, 5359, 5360, 5361, 5362, 5363, 5364, 5365, 5366, 5367, 5368, 5370, 5371, 5374, 5375, 5376, 5377, 5380, 5382, 5383, 5384,
5385, 5386, 5387, 5388, 5389, 5390, 5391, 5392, 5393, 5394, 5395, 5396, 5397, 5399, 5401, 5403, 5404, 5408, 5409, 5411, 5412, 5413, 5415, 5417, 5418,
5419, 5422, 5423, 5424, 5425, 5427, 5428, 5430, 5431, 5432, 5433, 5435, 5436, 5437, 5440, 5441, 5442, 5444, 5445, 5448, 5449, 5450, 5451, 5452, 5454,
5455, 5456, 5457, 5458, 5459, 5460, 5462, 5463, 5464, 5465, 5466, 5470, 5474, 5476, 5478, 5480, 5482, 5484, 5485, 5487, 5488, 5489, 5490, 5491, 5492,
5493, 5495, 5497, 5498, 5499, 5500, 5501, 5504, 5506, 5507, 5508, 5509, 5510, 5511, 5512, 5513, 5514, 5515, 5516, 5517, 5518, 5519, 5520, 5521, 5522,
5523, 5524, 5525, 5526, 5527, 5528, 5530, 5531, 5532, 5533, 5535, 5536, 5538, 5539, 5540, 5541, 5544, 5545, 5546, 5547, 5553, 5555, 5557, 5561, 5564,
5565, 5566, 5567, 5568, 5569, 5570, 5573, 5574, 5575, 5576, 5577, 5578, 5579, 5581, 5582, 5583, 5584, 5585, 5586, 5587, 5588, 5589, 5593, 5597, 5598,
5600, 5601, 5602, 5603, 5604, 5605, 5606, 5607, 5608, 5609, 5611, 5614, 5616, 5617, 5618, 5619, 5622, 5623, 5624, 5625, 5627, 5628, 5630, 5632, 5634,
5635, 5636, 5637, 5638, 5639, 5643, 5646, 5647, 5650, 5652, 5653, 5654, 5656, 5657, 5659, 5661, 5680, 5684, 5685, 5686, 5691, 5693, 5694, 5695, 5697,
5698, 5699, 5702, 5703, 5704, 5705, 5706, 5707, 5708, 5709, 5710, 5712, 5713, 5715, 5718, 5719, 5720, 5721, 5722, 5723, 5724, 5725, 5726, 5727, 5728,
5729, 5730, 5731, 5733, 5735, 5739, 5740, 5741, 5743, 5745, 5746, 5748, 5749, 5752, 5753, 5754, 5755, 5756, 5758, 5760, 5761, 5762, 5765, 5766, 5767,
5768, 5769, 5771, 5772, 5773, 5774, 5776, 5777, 5778, 5779, 5781, 5782, 5783, 5891, 5892, 5894, 5895, 5896, 5897, 5899, 5900, 5901, 5904, 5905, 5906,
5908, 5909, 5910, 5911, 5912, 5914, 5915, 5916, 5918, 5919, 5920, 5921, 5922, 5923, 5924, 5925, 5926, 5929, 5930, 5931, 5932, 5933, 5934, 5935, 5936,
5937, 5938, 5940, 5941, 5942, 5946, 5947, 5948, 5949, 5952, 5953, 5956, 5957, 5959, 5960, 5961, 5963, 5964, 5965, 5966, 5968, 5969, 5971, 5975, 5976,
5977, 5979, 5982, 5983, 5984, 5985, 5986, 5987, 5988, 5989, 5990, 5993, 5994, 5996, 5998, 6000, 6001, 6002, 6003, 6004, 6007, 6011, 6012, 6021, 6022,
6023, 6024, 6025, 6026, 6027, 6028, 6030, 6031, 6032, 6035, 6036, 6037, 6041, 6043, 6044, 6045, 6046, 6047, 6050, 6052, 6055, 6056, 6058, 6059, 6060,
6063, 6064, 6065, 6070, 6072, 6073, 6074, 6075, 6076, 6078, 6079, 6080, 6081, 6082, 6083, 6084, 6085, 6086, 6089, 6090, 6091, 6092, 6093, 6094, 6095,
6096, 6098, 6099, 6100, 6101, 6102, 6103, 6105, 6106, 6107, 6108, 6110, 6111, 6113, 6114, 6115, 6117, 6118, 6119, 6120, 6121, 6122, 6123, 6125, 6130,
6134, 6138, 6139, 6140, 6141, 6143, 6144, 6145, 6146, 6147, 6148, 6149, 6151, 6153, 6154, 6155, 6156, 6157, 6158, 6159, 6160, 6161, 6162, 6163, 6164,
6165, 6166, 6168, 6170, 6172, 6173, 6176, 6178, 6179, 6180, 6182, 6183, 6184, 6185, 6187, 6190, 6191, 6192, 6194, 6195, 6197, 6198, 6199, 6201, 6203,
6205, 6206, 6207, 6210, 6211, 6213, 6214, 6216, 6217, 6218, 6219, 6220, 6221, 6223, 6224, 6226, 6227, 6228, 6229, 6230, 6231, 6233, 6234, 6235, 6236,
6238, 6239, 6241, 6244, 6246, 6247, 6254, 6260, 6261, 6262, 6265, 6266, 6276, 6277, 6278, 6279, 6280, 6281, 6282, 6283, 6284, 6285, 6286, 6288, 6295,
6296, 6297, 6298, 6299, 6300, 6301, 6302, 6303, 6304, 6305, 6307, 6308, 6309, 6310, 6311, 6312, 6313, 6314, 6315, 6318, 6319, 6320, 6321, 6322, 6323,
6324, 6325, 6326, 6327, 6329, 6330, 6332, 6333, 6334, 6335, 6336, 6337, 6338, 6339, 6340, 6341, 6342, 6343, 6344, 6346, 6347, 6348, 6349, 6350, 6352,
6353, 6354, 6355, 6356, 6358, 6359, 6360, 6361, 6362, 6363, 6365, 6366, 6367, 6369, 6370, 6371, 6374, 6375, 6376, 6377, 6378, 6379, 6380, 6381, 6382,
6383, 6384, 6385, 6386, 6387, 6388, 6397, 6401, 6402, 6404, 6405, 6406, 6407, 6408, 6409, 6410, 6412, 6417, 6419, 6420, 6421, 6423, 6425, 6426, 6427,
6428, 6429, 6430, 6431, 6432, 6433, 6436, 6437, 6438, 6439, 6440, 6441, 6442, 6443, 6444, 6445, 6446, 6448, 6450, 6451, 6453, 6454, 6455, 6457, 6460,
6461, 6462, 6466, 6467, 6468, 6475, 6478, 6479, 6480, 6481, 6482, 6483, 6484, 6487, 6488, 6489, 6490, 6493, 6494, 6497, 6501, 6502, 6503, 6506, 6507,
6508, 6510, 6517, 6520, 6521, 6525, 6527, 6528, 6529, 6530, 6533, 6534, 6535, 6540, 6542, 6544, 6546, 6547, 6548, 6549, 6551, 6552, 6553, 6559, 6560,
6562, 6563, 6564, 6565, 6566, 6567, 6568, 6569, 6571, 6572, 6573, 6575, 6577, 6579, 6580, 6581, 6582, 6583, 6588, 6589, 6590, 6591, 6592, 6593, 6594,
6596, 6598, 6599, 6600, 6602, 6603, 6605, 6606, 6608, 6609, 6611, 6614, 6615, 6616, 6617, 6619, 6620, 6621, 6622, 6623, 6625, 6626, 6628, 6629, 6630,
6631, 6632, 6633, 6634, 6635, 6637, 6638, 6639, 6640, 6641, 6642, 6645, 6647, 6649, 6650, 6653, 6655, 6656, 6657, 6658, 6659, 6663, 6664, 6665, 6666,
6667, 6668, 6669, 6670, 6671, 6672, 6674, 6675, 6676, 6680, 6681, 6683, 6684, 6685, 6686, 6687, 6688, 6689, 6690, 6691, 6692, 6693, 6694, 6695, 6696,
6697, 6698, 6699, 6710, 6711, 6712, 6714, 6717, 6720, 6722, 6725, 6747, 6748, 6749, 6751, 6752, 6756, 6758, 6761, 6763, 6765, 6766, 6767, 6768, 6769,
6772, 6773, 6775, 6776, 6777, 6778, 6779, 6780, 6782, 6783, 6784, 6785, 6786, 6787, 6788, 6789, 6790, 6792, 6793, 6794, 6795, 6797, 6798, 6799, 6800,
6801, 6802, 6803, 6804, 6805, 6806, 6810, 6811, 6812, 6813, 6814, 6822, 6823, 6824, 6825, 6826, 6827, 6828, 6829, 6830, 6831, 6832, 6835, 6836, 6837,
6840, 6841, 6850, 6851, 6853, 6854, 6855, 6856, 6857, 6859, 6861, 6863, 6865, 6866, 6867, 6869, 6870, 6871, 6873, 6874, 6875, 6876, 6879, 6880, 6881,
6882, 6883, 6884, 6885, 6886, 6887, 6888, 6890, 6891, 6893, 6894, 6895, 6896, 6901, 6902, 6903, 6905, 6906, 6907, 6908, 6911, 6912, 6913, 6914, 6915,
6916, 6917, 6918, 6919, 6921, 6922, 6923, 6925, 6926, 6928, 6930, 6932, 6933, 6934, 6935, 6936, 6937, 6938, 6939, 6940, 6942, 6944, 6945, 6946, 6947,
6949, 6950, 6951, 6952, 6953, 6954, 6955, 6956, 6957, 6958, 6959, 6960, 6963, 6964, 6965, 6966, 6967, 6968, 6969, 6970, 6971, 6972, 6973, 6974, 6975,
6976, 6977, 6978, 6979, 6980, 6981, 6982, 6983, 6984, 6985, 6986, 6987, 6989, 6990, 6991, 6992, 6993, 6995, 6997, 6998, 6999, 7001, 7003, 7004, 7006,
7007, 7008, 7010, 7011, 7012, 7013, 7014, 7016, 7017, 7018, 7019, 7020, 7021, 7022, 7023, 7025, 7026, 7046, 7047, 7048, 7049, 7051, 7053, 7054, 7055,
7057, 7058, 7059, 7078, 7079, 7080, 7082, 7084, 7085, 7086, 7090, 7093, 7095, 7096, 7097, 7098, 7099, 7100, 7103, 7105, 7107, 7108, 7110, 7112, 7113,
7114, 7115, 7116, 7117, 7123, 7124, 7125, 7126, 7127, 7128, 7129, 7130, 7131, 7132, 7133, 7134, 7135, 7136, 7139, 7141, 7142, 7143, 7144, 7145, 7146,
7165, 7166, 7167, 7279, 7280, 7282, 7283, 7285, 7287, 7288, 7289, 7290, 7291, 7292, 7293, 7294, 7295, 7297, 7298, 7299, 7301, 7302, 7303, 7305, 7306,
7307, 7308, 7309, 7310, 7311, 7312, 7314, 7315, 7318, 7319, 7320, 7322, 7323, 7324, 7326, 7327, 7328, 7331, 7332, 7333, 7334, 7336, 7337, 7338, 7341,
7343, 7344, 7345, 7346, 7349, 7350, 7351, 7352, 7353, 7354, 7355, 7356, 7357, 7359, 7360, 7361, 7362, 7365, 7368, 7369, 7371, 7372, 7373, 7375, 7376,
7377, 7378, 7380, 7381, 7382, 7383, 7384, 7385, 7386, 7387, 7388, 7389, 7390, 7391, 7392, 7393, 7394, 7397, 7398, 7400, 7401, 7403, 7404, 7406, 7407,
7408, 7409, 7411, 7412, 7413, 7414, 7415, 7416, 7417, 7418, 7419, 7420, 7421, 7422, 7426, 7427, 7429, 7430, 7431, 7432, 7434, 7435, 7436, 7437, 7438,
7439, 7440, 7441, 7442, 7443, 7444, 7445, 7455, 7456, 7457, 7458, 7459, 7461, 7462, 7463, 7464, 7465, 7466, 7468, 7469, 7470, 7471, 7473, 7474, 7476,
7477, 7478, 7479, 7480, 7481, 7482, 7483, 7484, 7485, 7486, 7487, 7488, 7490, 7491, 7492, 7493, 7494, 7495, 7496, 7497, 7498, 7499, 7500, 7501, 7502,
7503, 7504, 7505, 7506, 7507, 7508, 7509, 7510, 7512, 7513, 7514, 7515, 7516, 7517, 7518, 7519, 7520, 7521, 7522, 7523, 7524, 7525, 7526, 7527, 7528,
7529, 7531, 7532, 7534, 7535, 7537, 7538, 7539, 7540, 7542, 7543, 7544, 7546, 7547, 7556, 7557, 7558, 7559, 7560, 7561, 7562, 7563, 7564, 7565, 7567,
7568, 7571, 7573, 7574, 7575, 7576, 7577, 7578, 7580, 7583, 7584, 7586, 7587, 7588, 7611, 7623, 7625, 7626, 7627, 7628, 7631, 7633, 7634, 7635, 7636,
7637, 7638, 7639, 7640, 7641, 7642, 7643, 7644, 7645, 7648, 7649, 7650, 7652, 7653, 7654, 7655, 7656, 7657, 7660, 7662, 7663, 7665, 7666, 7667, 7668,
7669, 7670, 7671, 7673, 7674, 7675, 7676, 7677, 7678, 7680, 7681, 7682, 7683, 7685, 7687, 7688, 7689, 7690, 7691, 7693, 7694, 7695, 7696, 7699, 7700,
7704, 7705, 7707, 7708, 7709, 7710, 7711, 7712, 7713, 7714, 7715, 7716, 7717, 7718, 7719, 7720, 7721, 7722, 7723, 7724, 7725, 7726, 7727, 7728, 7729,
7730, 7731, 7736, 7738, 7739, 7740, 7741, 7742, 7744, 7745, 7746, 7747, 7748, 7749, 7750, 7751, 7753, 7754, 7755, 7756, 7757, 7759, 7760, 7761, 7762,
7764, 7765, 7766, 7767, 7768, 7769, 7770};

  std::set<int> good_runlist_set(good_run_list_vec.begin(), good_run_list_vec.end());
  
  
  if (run_filter != 1) return 0;
  if (filter_level == 1){
    std::cout << "FHC mode" << std::endl;
    outfile_name = prefix_out + "_FHC.root";
  }else{
    std::cout << "RHC mode" << std::endl;
    outfile_name = prefix_out + "_RHC.root";
  }
  
  bool flag_data = true;

  TFile *file1 = new TFile(input_file);
  TTree *T_BDTvars = (TTree*)file1->Get("wcpselection/T_BDTvars");
  
  TTree *T_eval = (TTree*)file1->Get("wcpselection/T_eval");
  TTree *T_pot = (TTree*)file1->Get("wcpselection/T_pot");
  TTree *T_PFeval = (TTree*)file1->Get("wcpselection/T_PFeval");
  TTree *T_KINEvars = (TTree*)file1->Get("wcpselection/T_KINEvars");

  TFile *file2 = new TFile(outfile_name,"RECREATE");
  file2->mkdir("wcpselection");
  file2->cd("wcpselection");
  TTree *t4 = new TTree("T_BDTvars","T_BDTvars");
  TTree *t1 = new TTree("T_eval","T_eval");
  TTree *t2 = new TTree("T_pot","T_pot");
  TTree *t3 = new TTree("T_PFeval", "T_PFeval");
  TTree *t5 = new TTree("T_KINEvars", "T_KINEvars");

  EvalInfo eval;
  eval.file_type = new std::string();
  POTInfo pot;
  TaggerInfo tagger;
  PFevalInfo pfeval;
  KineInfo kine;

  kine.kine_energy_particle = new std::vector<float>;
  kine.kine_energy_info = new std::vector<int>;
  kine.kine_particle_type = new std::vector<int>;
  kine.kine_energy_included = new std::vector<int>;
    
  tagger.pio_2_v_dis2 = new std::vector<float>;
  tagger.pio_2_v_angle2 = new std::vector<float>;
  tagger.pio_2_v_acc_length = new std::vector<float>;
  tagger.pio_2_v_flag = new std::vector<float>;
  tagger.sig_1_v_angle = new std::vector<float>;
  tagger.sig_1_v_flag_single_shower = new std::vector<float>;
  tagger.sig_1_v_energy = new std::vector<float>;
  tagger.sig_1_v_energy_1 = new std::vector<float>;
  tagger.sig_1_v_flag = new std::vector<float>;
  tagger.sig_2_v_energy = new std::vector<float>;
  tagger.sig_2_v_shower_angle = new std::vector<float>;
  tagger.sig_2_v_flag_single_shower = new std::vector<float>;
  tagger.sig_2_v_medium_dQ_dx = new std::vector<float>;
  tagger.sig_2_v_start_dQ_dx = new std::vector<float>;
  tagger.sig_2_v_flag = new std::vector<float>;
  tagger.stw_2_v_medium_dQ_dx = new std::vector<float>;
  tagger.stw_2_v_energy = new std::vector<float>;
  tagger.stw_2_v_angle = new std::vector<float>;
  tagger.stw_2_v_dir_length = new std::vector<float>;
  tagger.stw_2_v_max_dQ_dx = new std::vector<float>;
  tagger.stw_2_v_flag = new std::vector<float>;
  tagger.stw_3_v_angle = new std::vector<float>;
  tagger.stw_3_v_dir_length = new std::vector<float>;
  tagger.stw_3_v_energy = new std::vector<float>;
  tagger.stw_3_v_medium_dQ_dx = new std::vector<float>;
  tagger.stw_3_v_flag = new std::vector<float>;
  tagger.stw_4_v_angle = new std::vector<float>;
  tagger.stw_4_v_dis = new std::vector<float>;
  tagger.stw_4_v_energy = new std::vector<float>;
  tagger.stw_4_v_flag = new std::vector<float>;
  tagger.br3_3_v_energy = new std::vector<float>;
  tagger.br3_3_v_angle = new std::vector<float>;
  tagger.br3_3_v_dir_length = new std::vector<float>;
  tagger.br3_3_v_length = new std::vector<float>;
  tagger.br3_3_v_flag = new std::vector<float>;
  tagger.br3_5_v_dir_length = new std::vector<float>;
  tagger.br3_5_v_total_length = new std::vector<float>;
  tagger.br3_5_v_flag_avoid_muon_check = new std::vector<float>;
  tagger.br3_5_v_n_seg = new std::vector<float>;
  tagger.br3_5_v_angle = new std::vector<float>;
  tagger.br3_5_v_sg_length = new std::vector<float>;
  tagger.br3_5_v_energy = new std::vector<float>;
  tagger.br3_5_v_n_main_segs = new std::vector<float>;
  tagger.br3_5_v_n_segs = new std::vector<float>;
  tagger.br3_5_v_shower_main_length = new std::vector<float>;
  tagger.br3_5_v_shower_total_length = new std::vector<float>;
  tagger.br3_5_v_flag = new std::vector<float>;
  tagger.br3_6_v_angle = new std::vector<float>;
  tagger.br3_6_v_angle1 = new std::vector<float>;
  tagger.br3_6_v_flag_shower_trajectory = new std::vector<float>;
  tagger.br3_6_v_direct_length = new std::vector<float>;
  tagger.br3_6_v_length = new std::vector<float>;
  tagger.br3_6_v_n_other_vtx_segs = new std::vector<float>;
  tagger.br3_6_v_energy = new std::vector<float>;
  tagger.br3_6_v_flag = new std::vector<float>;
  tagger.tro_1_v_particle_type = new std::vector<float>;
  tagger.tro_1_v_flag_dir_weak = new std::vector<float>;
  tagger.tro_1_v_min_dis = new std::vector<float>;
  tagger.tro_1_v_sg1_length = new std::vector<float>;
  tagger.tro_1_v_shower_main_length = new std::vector<float>;
  tagger.tro_1_v_max_n_vtx_segs = new std::vector<float>;
  tagger.tro_1_v_tmp_length = new std::vector<float>;
  tagger.tro_1_v_medium_dQ_dx = new std::vector<float>;
  tagger.tro_1_v_dQ_dx_cut = new std::vector<float>;
  tagger.tro_1_v_flag_shower_topology = new std::vector<float>;
  tagger.tro_1_v_flag = new std::vector<float>;
  tagger.tro_2_v_energy = new std::vector<float>;
  tagger.tro_2_v_stem_length = new std::vector<float>;
  tagger.tro_2_v_iso_angle = new std::vector<float>;
  tagger.tro_2_v_max_length = new std::vector<float>;
  tagger.tro_2_v_angle = new std::vector<float>;
  tagger.tro_2_v_flag = new std::vector<float>;
  tagger.tro_4_v_dir2_mag = new std::vector<float>;
  tagger.tro_4_v_angle = new std::vector<float>;
  tagger.tro_4_v_angle1 = new std::vector<float>;
  tagger.tro_4_v_angle2 = new std::vector<float>;
  tagger.tro_4_v_length = new std::vector<float>;
  tagger.tro_4_v_length1 = new std::vector<float>;
  tagger.tro_4_v_medium_dQ_dx = new std::vector<float>;
  tagger.tro_4_v_end_dQ_dx = new std::vector<float>;
  tagger.tro_4_v_energy = new std::vector<float>;
  tagger.tro_4_v_shower_main_length = new std::vector<float>;
  tagger.tro_4_v_flag_shower_trajectory = new std::vector<float>;
  tagger.tro_4_v_flag = new std::vector<float>;
  tagger.tro_5_v_max_angle = new std::vector<float>;
  tagger.tro_5_v_min_angle = new std::vector<float>;
  tagger.tro_5_v_max_length = new std::vector<float>;
  tagger.tro_5_v_iso_angle = new std::vector<float>;
  tagger.tro_5_v_n_vtx_segs = new std::vector<float>;
  tagger.tro_5_v_min_count = new std::vector<float>;
  tagger.tro_5_v_max_count = new std::vector<float>;
  tagger.tro_5_v_energy = new std::vector<float>;
  tagger.tro_5_v_flag = new std::vector<float>;
  tagger.lol_1_v_energy = new std::vector<float>;
  tagger.lol_1_v_vtx_n_segs = new std::vector<float>;
  tagger.lol_1_v_nseg = new std::vector<float>;
  tagger.lol_1_v_angle = new std::vector<float>;
  tagger.lol_1_v_flag = new std::vector<float>;
  tagger.lol_2_v_length = new std::vector<float>;
  tagger.lol_2_v_angle = new std::vector<float>;
  tagger.lol_2_v_type = new std::vector<float>;
  tagger.lol_2_v_vtx_n_segs = new std::vector<float>;
  tagger.lol_2_v_energy = new std::vector<float>;
  tagger.lol_2_v_shower_main_length = new std::vector<float>;
  tagger.lol_2_v_flag_dir_weak = new std::vector<float>;
  tagger.lol_2_v_flag = new std::vector<float>;
  tagger.cosmict_flag_10 = new std::vector<float>;
  tagger.cosmict_10_flag_inside = new std::vector<float>;
  tagger.cosmict_10_vtx_z = new std::vector<float>;
  tagger.cosmict_10_flag_shower = new std::vector<float>;
  tagger.cosmict_10_flag_dir_weak = new std::vector<float>;
  tagger.cosmict_10_angle_beam = new std::vector<float>;
  tagger.cosmict_10_length = new std::vector<float>;
  tagger.numu_cc_flag_1 = new std::vector<float>;
  tagger.numu_cc_1_particle_type = new std::vector<float>;
  tagger.numu_cc_1_length = new std::vector<float>;
  tagger.numu_cc_1_medium_dQ_dx = new std::vector<float>;
  tagger.numu_cc_1_dQ_dx_cut = new std::vector<float>;
  tagger.numu_cc_1_direct_length = new std::vector<float>;
  tagger.numu_cc_1_n_daughter_tracks = new std::vector<float>;
  tagger.numu_cc_1_n_daughter_all = new std::vector<float>;
  tagger.numu_cc_flag_2 = new std::vector<float>;
  tagger.numu_cc_2_length = new std::vector<float>;
  tagger.numu_cc_2_total_length = new std::vector<float>;
  tagger.numu_cc_2_n_daughter_tracks = new std::vector<float>;
  tagger.numu_cc_2_n_daughter_all = new std::vector<float>;
  tagger.pio_2_v_dis2 = new std::vector<float>;
  tagger.pio_2_v_angle2 = new std::vector<float>;
  tagger.pio_2_v_acc_length = new std::vector<float>;
  tagger.pio_2_v_flag = new std::vector<float>;
  tagger.sig_1_v_angle = new std::vector<float>;
  tagger.sig_1_v_flag_single_shower = new std::vector<float>;
  tagger.sig_1_v_energy = new std::vector<float>;
  tagger.sig_1_v_energy_1 = new std::vector<float>;
  tagger.sig_1_v_flag = new std::vector<float>;
  tagger.sig_2_v_energy = new std::vector<float>;
  tagger.sig_2_v_shower_angle = new std::vector<float>;
  tagger.sig_2_v_flag_single_shower = new std::vector<float>;
  tagger.sig_2_v_medium_dQ_dx = new std::vector<float>;
  tagger.sig_2_v_start_dQ_dx = new std::vector<float>;
  tagger.sig_2_v_flag = new std::vector<float>;
  tagger.stw_2_v_medium_dQ_dx = new std::vector<float>;
  tagger.stw_2_v_energy = new std::vector<float>;
  tagger.stw_2_v_angle = new std::vector<float>;
  tagger.stw_2_v_dir_length = new std::vector<float>;
  tagger.stw_2_v_max_dQ_dx = new std::vector<float>;
  tagger.stw_2_v_flag = new std::vector<float>;
  tagger.stw_3_v_angle = new std::vector<float>;
  tagger.stw_3_v_dir_length = new std::vector<float>;
  tagger.stw_3_v_energy = new std::vector<float>;
  tagger.stw_3_v_medium_dQ_dx = new std::vector<float>;
  tagger.stw_3_v_flag = new std::vector<float>;
  tagger.stw_4_v_angle = new std::vector<float>;
  tagger.stw_4_v_dis = new std::vector<float>;
  tagger.stw_4_v_energy = new std::vector<float>;
  tagger.stw_4_v_flag = new std::vector<float>;
  tagger.br3_3_v_energy = new std::vector<float>;
  tagger.br3_3_v_angle = new std::vector<float>;
  tagger.br3_3_v_dir_length = new std::vector<float>;
  tagger.br3_3_v_length = new std::vector<float>;
  tagger.br3_3_v_flag = new std::vector<float>;
  tagger.br3_5_v_dir_length = new std::vector<float>;
  tagger.br3_5_v_total_length = new std::vector<float>;
  tagger.br3_5_v_flag_avoid_muon_check = new std::vector<float>;
  tagger.br3_5_v_n_seg = new std::vector<float>;
  tagger.br3_5_v_angle = new std::vector<float>;
  tagger.br3_5_v_sg_length = new std::vector<float>;
  tagger.br3_5_v_energy = new std::vector<float>;
  tagger.br3_5_v_n_main_segs = new std::vector<float>;
  tagger.br3_5_v_n_segs = new std::vector<float>;
  tagger.br3_5_v_shower_main_length = new std::vector<float>;
  tagger.br3_5_v_shower_total_length = new std::vector<float>;
  tagger.br3_5_v_flag = new std::vector<float>;
  tagger.br3_6_v_angle = new std::vector<float>;
  tagger.br3_6_v_angle1 = new std::vector<float>;
  tagger.br3_6_v_flag_shower_trajectory = new std::vector<float>;
  tagger.br3_6_v_direct_length = new std::vector<float>;
  tagger.br3_6_v_length = new std::vector<float>;
  tagger.br3_6_v_n_other_vtx_segs = new std::vector<float>;
  tagger.br3_6_v_energy = new std::vector<float>;
  tagger.br3_6_v_flag = new std::vector<float>;
  tagger.tro_1_v_particle_type = new std::vector<float>;
  tagger.tro_1_v_flag_dir_weak = new std::vector<float>;
  tagger.tro_1_v_min_dis = new std::vector<float>;
  tagger.tro_1_v_sg1_length = new std::vector<float>;
  tagger.tro_1_v_shower_main_length = new std::vector<float>;
  tagger.tro_1_v_max_n_vtx_segs = new std::vector<float>;
  tagger.tro_1_v_tmp_length = new std::vector<float>;
  tagger.tro_1_v_medium_dQ_dx = new std::vector<float>;
  tagger.tro_1_v_dQ_dx_cut = new std::vector<float>;
  tagger.tro_1_v_flag_shower_topology = new std::vector<float>;
  tagger.tro_1_v_flag = new std::vector<float>;
  tagger.tro_2_v_energy = new std::vector<float>;
  tagger.tro_2_v_stem_length = new std::vector<float>;
  tagger.tro_2_v_iso_angle = new std::vector<float>;
  tagger.tro_2_v_max_length = new std::vector<float>;
  tagger.tro_2_v_angle = new std::vector<float>;
  tagger.tro_2_v_flag = new std::vector<float>;
  tagger.tro_4_v_dir2_mag = new std::vector<float>;
  tagger.tro_4_v_angle = new std::vector<float>;
  tagger.tro_4_v_angle1 = new std::vector<float>;
  tagger.tro_4_v_angle2 = new std::vector<float>;
  tagger.tro_4_v_length = new std::vector<float>;
  tagger.tro_4_v_length1 = new std::vector<float>;
  tagger.tro_4_v_medium_dQ_dx = new std::vector<float>;
  tagger.tro_4_v_end_dQ_dx = new std::vector<float>;
  tagger.tro_4_v_energy = new std::vector<float>;
  tagger.tro_4_v_shower_main_length = new std::vector<float>;
  tagger.tro_4_v_flag_shower_trajectory = new std::vector<float>;
  tagger.tro_4_v_flag = new std::vector<float>;
  tagger.tro_5_v_max_angle = new std::vector<float>;
  tagger.tro_5_v_min_angle = new std::vector<float>;
  tagger.tro_5_v_max_length = new std::vector<float>;
  tagger.tro_5_v_iso_angle = new std::vector<float>;
  tagger.tro_5_v_n_vtx_segs = new std::vector<float>;
  tagger.tro_5_v_min_count = new std::vector<float>;
  tagger.tro_5_v_max_count = new std::vector<float>;
  tagger.tro_5_v_energy = new std::vector<float>;
  tagger.tro_5_v_flag = new std::vector<float>;
  tagger.lol_1_v_energy = new std::vector<float>;
  tagger.lol_1_v_vtx_n_segs = new std::vector<float>;
  tagger.lol_1_v_nseg = new std::vector<float>;
  tagger.lol_1_v_angle = new std::vector<float>;
  tagger.lol_1_v_flag = new std::vector<float>;
  tagger.lol_2_v_length = new std::vector<float>;
  tagger.lol_2_v_angle = new std::vector<float>;
  tagger.lol_2_v_type = new std::vector<float>;
  tagger.lol_2_v_vtx_n_segs = new std::vector<float>;
  tagger.lol_2_v_energy = new std::vector<float>;
  tagger.lol_2_v_shower_main_length = new std::vector<float>;
  tagger.lol_2_v_flag_dir_weak = new std::vector<float>;
  tagger.lol_2_v_flag = new std::vector<float>;
  tagger.cosmict_flag_10 = new std::vector<float>;
  tagger.cosmict_10_flag_inside = new std::vector<float>;
  tagger.cosmict_10_vtx_z = new std::vector<float>;
  tagger.cosmict_10_flag_shower = new std::vector<float>;
  tagger.cosmict_10_flag_dir_weak = new std::vector<float>;
  tagger.cosmict_10_angle_beam = new std::vector<float>;
  tagger.cosmict_10_length = new std::vector<float>;
  tagger.numu_cc_flag_1 = new std::vector<float>;
  tagger.numu_cc_1_particle_type = new std::vector<float>;
  tagger.numu_cc_1_length = new std::vector<float>;
  tagger.numu_cc_1_medium_dQ_dx = new std::vector<float>;
  tagger.numu_cc_1_dQ_dx_cut = new std::vector<float>;
  tagger.numu_cc_1_direct_length = new std::vector<float>;
  tagger.numu_cc_1_n_daughter_tracks = new std::vector<float>;
  tagger.numu_cc_1_n_daughter_all = new std::vector<float>;
  tagger.numu_cc_flag_2 = new std::vector<float>;
  tagger.numu_cc_2_length = new std::vector<float>;
  tagger.numu_cc_2_total_length = new std::vector<float>;
  tagger.numu_cc_2_n_daughter_tracks = new std::vector<float>;
  tagger.numu_cc_2_n_daughter_all = new std::vector<float>;

  set_tree_address(T_BDTvars, tagger,2 );
  put_tree_address(t4, tagger,2);

  if (flag_data){
    set_tree_address(T_eval, eval,2);
    put_tree_address(t1, eval,2);

    set_tree_address(T_PFeval, pfeval,2);
    put_tree_address(t3, pfeval,2);
  }else{
    set_tree_address(T_eval, eval);
    put_tree_address(t1, eval);

    set_tree_address(T_PFeval, pfeval);
    put_tree_address(t3, pfeval);
  }

  set_tree_address(T_pot, pot);
  put_tree_address(t2, pot);

  set_tree_address(T_KINEvars, kine);
  put_tree_address(t5, kine);
    

  T_eval->SetBranchStatus("*",1);
  T_BDTvars->SetBranchStatus("*",1);


  
  for (int i=0;i!=T_BDTvars->GetEntries();i++){
    T_BDTvars->GetEntry(i);
    T_eval->GetEntry(i); tagger.match_isFC = eval.match_isFC;
    T_KINEvars->GetEntry(i); tagger.kine_reco_Enu = kine.kine_reco_Enu;
    T_PFeval->GetEntry(i);

    if (good_runlist_set.find(eval.run) == good_runlist_set.end()) continue;
    
    if (filter_level==1 && eval.run > 6748 && eval.run <=7001) continue;
    if (filter_level!=1 && eval.run <=6748) continue;
    // if (flag_remove){
    //   clear_tagger_info(tagger);
    //   clear_kine_info(kine);
    //   clear_pfeval_info(pfeval);
    // }
    
    t4->Fill();
    t1->Fill();
    t3->Fill();
    t5->Fill();
  }

  for (Int_t i=0;i!=T_pot->GetEntries();i++){
    T_pot->GetEntry(i);

    if (good_runlist_set.find(pot.runNo) == good_runlist_set.end()) continue;
    
    if (filter_level==1 && pot.runNo > 6748 && pot.runNo <=7001) continue;
    if (filter_level!=1 && pot.runNo <=6748) continue;
    t2->Fill();
  }


  
  file2->Write();
  file2->Close();


  return 0;

  
  
}
