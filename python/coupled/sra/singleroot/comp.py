"""
    compare to values from the Matlab implementation 
"""

import numpy as np
import matplotlib.pyplot as plt

dx = 1.
z_ = np.linspace(0 - dx / 2, -100 + dx / 2, 100)  # cell centers

jan_kxup = np.loadtxt("Kxlinupscale_singlerootconstk100.txt")
jan_suf = np.loadtxt("SUFupscale_singlerootconstk100.txt")
print("sum suf", np.sum(jan_suf))
print(list(jan_kxup))
print(list(jan_suf))

suf = np.array([0.023912050833266352, 0.02336006009061267, 0.022821279773238265, 0.022295405193693822, 0.021782138962896783, 0.02128119082195379, 0.020792277478015476, 0.02031512244407053, 0.01984945588258866, 0.019395014452923836, 0.018951541162391658, 0.01851878522093657, 0.018096501899306745, 0.017684452390656406, 0.017282403675497424, 0.016890128389923636, 0.016507404697033556, 0.016134016161478647, 0.01576975162706626, 0.015414405097347945, 0.015067775619125784, 0.014729667168810622, 0.014399888541568153, 0.01407825324319001, 0.013764579384628823, 0.013458689579137536, 0.013160410841954831, 0.012869574492479942, 0.01258601605888153, 0.012309575185086684, 0.012040095540097384, 0.01177742472958329, 0.011521414209700683, 0.011271919203088955, 0.011028798616997074, 0.010791914963493762, 0.010561134281716232, 0.010336326062113533, 0.010117363172641676, 0.009904121786868737, 0.009696481313949398, 0.009494324330429184, 0.009297536513839959, 0.009106006578049027, 0.008919626210325358, 0.008738290010087281, 0.008561895429297065, 0.008390342714468622, 0.008223534850255591, 0.008061377504587868, 0.00790377897532555, 0.007750650138400174, 0.007601904397413848, 0.0074574576346678485, 0.007317228163592934, 0.0071811366825544985, 0.0070491062300064525, 0.006921062140968433, 0.006796932004801757, 0.006676645624260261, 0.0065601349757928, 0.006447334171075033, 0.006338179419748702, 0.006232608993347328, 0.006130563190387942, 0.006031984302609094, 0.00593681658233608, 0.005845006210954872, 0.005756501268476981, 0.005671251704178011, 0.00558920930829331, 0.0055103276847547015, 0.005434562224952889, 0.005361870082510701, 0.005292210149052885, 0.005225543030958777, 0.005161831027084681, 0.005101038107443378, 0.005043129892828677, 0.004988073635373523, 0.004935838200030643, 0.004886394046965256, 0.004839713214849907, 0.004795769305051962, 0.0047545374667048316, 0.004715994382654471, 0.004680118256273221, 0.004646888799133525, 0.004616287219534553, 0.004588296211875241, 0.004562899946867747, 0.004540084062585778, 0.0045198356563427265, 0.004502143277395031, 0.004486996920466623, 0.004474388020090821, 0.004464309445766435, 0.004456755497925375, 0.0044517219047094745, 0.004449205819554691])
# suf = np.array([0.023633827680402498, 0.023088493548372863, 0.02255621626728385, 0.022036694827239824, 0.021529635431943264, 0.021034751332549013, 0.020551762665504004, 0.020080396294280663, 0.01962038565491461, 0.019171470605259134, 0.018733397277871384, 0.01830591793644692, 0.017888790835721546, 0.017481780084761134, 0.017084655513562152, 0.01669719254288748, 0.016319172057263832, 0.015950380281069064, 0.015590608657639183, 0.015239653731326793, 0.014897317032444214, 0.01456340496502622, 0.014237728697348946, 0.013920104055143028, 0.013610351417440618, 0.013308295614997304, 0.01301376583123158, 0.012726595505625804, 0.012446622239533984, 0.012173687704343206, 0.011907637551936667, 0.011648321327407776, 0.011395592383975897, 0.01114930780005563, 0.010909328298432762, 0.010675518167501138, 0.010447745184515946, 0.01022588054081998, 0.010009798769000652, 0.00979937767193647, 0.009594498253692955, 0.009395044652228856, 0.009200904073874607, 0.009011966729546015, 0.008828125772657068, 0.008649277238696768, 0.008475319986435808, 0.00830615564072987, 0.008141688536887157, 0.007981825666568743, 0.007826476625191117, 0.007675553560801189, 0.007528971124394832, 0.007386646421650888, 0.007248498966053326, 0.007114450633375043, 0.006984425617497568, 0.0068583503875417, 0.0067361536462848095, 0.006617766289841313, 0.006503121368583502, 0.006392154049280636, 0.006284801578434891, 0.00618100324679342, 0.00608070035501645, 0.0059838361804820375, 0.005890355945208664, 0.0058002067848775575, 0.005713337718937215, 0.005629699621773225, 0.005549245194927072, 0.0054719289403482325, 0.00539770713466442, 0.005326537804455434, 0.005258380702516633, 0.005193197285098605, 0.0051309506901101596, 0.005071605716272328, 0.005015128803211567, 0.004961488012480927, 0.004910653009498426, 0.004862595046392448, 0.004817286945744433, 0.004774703085219684, 0.004734819383077599, 0.004697613284553121, 0.004663063749101726, 0.004631151238500705, 0.004601857705800045, 0.004575166585116633, 0.004551062782266038, 0.0045295326662265385, 0.0045105640614306, 0.0044941462408794345, 0.0044802699200767356, 0.004468927251778171, 0.004460111821553666, 0.004453818644159954, 0.004450044160721355, 0.004448786236717184])
kxup = np.array([0.0048068120307766, 0.0015852421579424256, 0.0009405529297225264, 0.0006642814887913048, 0.0005108414047982881, 0.00041323954866687537, 0.0003457064012061989, 0.0002962153407210578, 0.00025839902548839105, 0.00022857097042964947, 0.00020444900503736074, 0.0001845447264821264, 0.00016784600439320732, 0.00015364058979695208, 0.0001414127111364977, 0.00013077969644833566, 0.00012145164137781665, 0.00011320490494164761, 0.00010586420104774593, 9.92902001892111e-05, 9.337075983207251e-05, 8.801460208712262e-05, 8.314667730430229e-05, 7.870471140991043e-05, 7.463659878087657e-05, 7.08984085521812e-05, 6.745284232139807e-05, 6.426802835129684e-05, 6.131656962329384e-05, 5.8574785506993874e-05, 5.6022102610112864e-05, 5.364056165816905e-05, 5.141441541336922e-05, 4.932979860993718e-05, 4.7374455294150957e-05, 4.553751225010147e-05, 4.380928967303324e-05, 4.218114213757966e-05, 4.06453243529073e-05, 3.919487731231742e-05, 3.782353131249431e-05, 3.6525622997068444e-05, 3.529602411475627e-05, 3.413008010711812e-05, 3.302355697984903e-05, 3.197259518335903e-05, 3.0973669447630424e-05, 3.0023553694029943e-05, 2.9119290291474002e-05, 2.825816304275148e-05, 2.743767338411426e-05, 2.6655519361547338e-05, 2.5909577013665617e-05, 2.5197883846529113e-05, 2.4518624131872237e-05, 2.3870115798956768e-05, 2.3250798722804238e-05, 2.2659224239017272e-05, 2.2094045738629744e-05, 2.155401021614294e-05, 2.103795066068995e-05, 2.054477919459987e-05, 2.0073480875900527e-05, 1.9623108091827256e-05, 1.9192775479466434e-05, 1.8781655317479425e-05, 1.8388973339612614e-05, 1.801400492655884e-05, 1.76560716378271e-05, 1.7314538049710177e-05, 1.6988808869308512e-05, 1.667832629795158e-05, 1.6382567620322672e-05, 1.610104299819613e-05, 1.5833293449986386e-05, 1.5578888999327497e-05, 1.533742697768541e-05, 1.5108530467583767e-05, 1.489184687442376e-05, 1.4687046616122482e-05, 1.4493821920901153e-05, 1.4311885724542447e-05, 1.4140970659318991e-05, 1.3980828127585767e-05, 1.3831227453739318e-05, 1.369195510888531e-05, 1.3562814003132377e-05, 1.3443622840951349e-05, 1.3334215535511557e-05, 1.3234440678336015e-05, 1.314416106100979e-05, 1.3063253246035313e-05, 1.2991607184259267e-05, 1.2929125876601065e-05, 1.28757250780966e-05, 1.2831333042535634e-05, 1.279589030621957e-05, 1.276934950960108e-05, 1.2751675255790379e-05, 1.2742844005126687e-05])
# kxup = np.array([0.0023885381408769125, 0.001181703854401832, 0.0007791160135063994, 0.0005778117213528334, 0.0004570565420917627, 0.00037658593829532444, 0.0003191388472069964, 0.00027608315206986624, 0.0002426226066938632, 0.00021587918235454014, 0.00019402126745173456, 0.00017582771392080437, 0.00016045307023370144, 0.00014729341826239434, 0.00013590587013881538, 0.00012595825102088994, 0.00011719653934931576, 0.0001094231594685911, 0.00010248212868241161, 9.624865986792274e-05, 9.062173456664762e-05, 8.551870146369686e-05, 8.087128387224917e-05, 7.662258528870619e-05, 7.272481357526425e-05, 6.913753030432682e-05, 6.582628911919617e-05, 6.27615658621763e-05, 5.991791004713469e-05, 5.727326603211105e-05, 5.480842557452012e-05, 5.250658303042598e-05, 5.035297142597182e-05, 4.8334562751636884e-05, 4.6439819635398303e-05, 4.46584884049679e-05, 4.29814257091575e-05, 4.140045251677655e-05, 3.9908230579444466e-05, 3.849815742739158e-05, 3.716427673431214e-05, 3.590120148997149e-05, 3.470404789578021e-05, 3.356837827759395e-05, 3.2490151613237416e-05, 3.146568051616298e-05, 3.0491593713862328e-05, 2.9564803219882597e-05, 2.8682475529101116e-05, 2.7842006273174016e-05, 2.7040997861406096e-05, 2.6277239705335e-05, 2.5548690685965e-05, 2.4853463573122086e-05, 2.4189811148668885e-05, 2.3556113820793664e-05, 2.29508685464618e-05, 2.2372678904358297e-05, 2.1820246182039612e-05, 2.1292361359196484e-05, 2.0787897884429462e-05, 2.0305805156190963e-05, 1.9845102629906805e-05, 1.940487448305303e-05, 1.8984264778376182e-05, 1.8582473072711222e-05, 1.8198750425142197e-05, 1.7832395763710762e-05, 1.7482752574626203e-05, 1.7149205882070345e-05, 1.683117949030658e-05, 1.6528133462968387e-05, 1.6239561817179135e-05, 1.596499041259614e-05, 1.5703975017621652e-05, 1.545609953692027e-05, 1.5220974386059937e-05, 1.4998235000579329e-05, 1.4787540468103649e-05, 1.458857227330369e-05, 1.440103314653851e-05, 1.4224646007955398e-05, 1.4059152999656176e-05, 1.3904314599287642e-05, 1.3759908809087537e-05, 1.3625730415023423e-05, 1.3501590311210094e-05, 1.3387314885286886e-05, 1.328274546088713e-05, 1.318773779374242e-05, 1.3102161618339987e-05, 1.3025900242395888e-05, 1.2958850186724572e-05, 1.2900920868379455e-05, 1.2852034325212769e-05, 1.2812124980259354e-05, 1.2781139444590205e-05, 1.2759036357510248e-05, 1.2745786263193148e-05, 1.2741371523055983e-05])
print("sum suf", np.sum(suf))

plt.plot(jan_suf, z_, label = "matlab")
plt.plot(suf, z_, label = "cplantbox")
plt.xlabel("SUF (1)")
plt.ylabel("depth (cm)")
plt.legend()
plt.show()

plt.plot(jan_kxup, z_, label = "matlab")
plt.plot(kxup, z_, label = "cplantbox")
plt.xlabel("kx_up (cm2/day)")
plt.ylabel("depth (cm)")
plt.legend()
plt.show()

krs = 0.002337363248476863
kr_surf = 1.8e-4 * 2 * 0.05 * np.pi * 1
print("krs", krs)
print("kr_surf", kr_surf)
print()

print("matlab")
print("suf", jan_suf[0])
suf_krs = krs * jan_suf[0]
print("suf_krs", suf_krs)
print("1. / (kr_surf - suf_krs)", 1. / (kr_surf - suf_krs))
print()

print("cplantbox")
print("suf", suf[0])
suf_krs = krs * suf[0]
print("suf_krs", suf_krs)
print("1. / (kr_surf - suf_krs)", 1. / (kr_surf - suf_krs))
print()
