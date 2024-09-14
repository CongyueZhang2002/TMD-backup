      PROGRAM EIC
      IMPLICIT NONE
      integer nloops,hop,nll,prepdf,preff
      common /scheme/ nloops,hop,nll,prepdf,preff
      REAL*8 aN, bN, g2A, A2, g2B, B2
      COMMON /FITP/ aN, bN, g2A, A2, g2B, B2
      real*8 Ecm,xb,zh,Q0,Q2,Q,PhT,APFEL_init
      double precision eps
      integer IH,IC,IT
      integer H_AA
      COMMON /HMASS/ H_AA
!      COMMON /REPLICA_INFO/ REPLICA_NUM
      COMMON /meson/ IH,IC
      real*8 F2light_D,FLlight_D,disD
      real*8 F2light_A,FLlight_A,disA
      real*8 fuuA,fuuD
      integer i,j,k
      real*8 fu,fub,fd,fdb,fs,fsb,fg,fC,fCB,fB,fBB
      real*8 F2light,FLlight
      real*8 aNs(201),bNs(201)
CC----Weird Parameter which APFEL uses --> I don't care to understand precisely why its needed (I have an idea) but I will just use it.
      parameter(eps=1d-10)
      data aNs /
     >  0.01706923088211247,
     >  0.021866956197765325,
     >  0.017815406077917054,
     >  0.01987487028429848,
     >  0.017585099841989527,
     >  0.012888401818107714,
     >  0.01781576688748901,
     >  0.013963863856388205,
     >  0.01376816186851837,
     >  0.018205346385391688,
     >  0.018541737255177448,
     >  0.016117366272188738,
     >  0.01362824446958276,
     >  0.017449149002092046,
     >  0.019274505121635,
     >  0.02240817720563477,
     >  0.0183247913971832,
     >  0.013622087845315434,
     >  0.018321965674514287,
     >  0.021206353078312826,
     >  0.02156375107030807,
     >  0.018859982474512427,
     >  0.017431980256055454,
     >  0.018324742953981146,
     >  0.01746461745291028,
     >  0.015972460007468894,
     >  0.015977379900557958,
     >  0.01566337956886626,
     >  0.01862419263065029,
     >  0.015977730955473043,
     >  0.014854577629987915,
     >  0.01742855981000627,
     >  0.019131172742019725,
     >  0.018109551909356408,
     >  0.019875417585738398,
     >  0.015825570356174212,
     >  0.017911793604347137,
     >  0.015563836066546809,
     >  0.01393723038359708,
     >  0.016087633808475527,
     >  0.01979198156137996,
     >  0.011031625005990949,
     >  0.015820521406463086,
     >  0.022689362518688444,
     >  0.019884453244182032,
     >  0.017459047681419038,
     >  0.01832469635056037,
     >  0.015290370366882687,
     >  0.014853018579429375,
     >  0.01582863148072863,
     >  0.01682400951636684,
     >  0.01832390233994094,
     >  0.022643085100038027,
     >  0.015977524997301554,
     >  0.01825086193174716,
     >  0.0204111389756005,
     >  0.019358178766246754,
     >  0.014855921300641615,
     >  0.018842622587408226,
     >  0.015405200998371588,
     >  0.014547140018816157,
     >  0.013923518577912145,
     >  0.021028858033007608,
     >  0.013926376614824584,
     >  0.019488812177874775,
     >  0.018822312707757846,
     >  0.015562784928725525,
     >  0.01725213043378865,
     >  0.018820585167140687,
     >  0.01949874470913586,
     >  0.020466961335998485,
     >  0.01988371459438437,
     >  0.018205552977482387,
     >  0.018254169526112737,
     >  0.01832428346449061,
     >  0.015703175009508474,
     >  0.017432193796188822,
     >  0.013213958039744178,
     >  0.019201858778184425,
     >  0.017440187633876225,
     >  0.018321057880951702,
     >  0.01944756450126138,
     >  0.01927688772012999,
     >  0.012246993045045285,
     >  0.015766916710201282,
     >  0.014486213807759147,
     >  0.011601447097327078,
     >  0.015586045121661438,
     >  0.017033904806391942,
     >  0.020373735304603846,
     >  0.021145753915987872,
     >  0.019662620596988426,
     >  0.017125868048177586,
     >  0.021506026852679508,
     >  0.01884731916468833,
     >  0.018373988381781527,
     >  0.017427923973416155,
     >  0.01687988762430287,
     >  0.015842294028822715,
     >  0.017570300972535813,
     >  0.021455982893284736,
     >  0.012723335945112088,
     >  0.019662409727333297,
     >  0.013715621424243984,
     >  0.017818329909962183,
     >  0.01448608510372941,
     >  0.014120798959328731,
     >  0.017444467681592857,
     >  0.015135774417141604,
     >  0.01782095706544562,
     >  0.01570730629401726,
     >  0.015977720333270037,
     >  0.013487868054093593,
     >  0.015386564762736116,
     >  0.01782210208134774,
     >  0.015630050546936242,
     >  0.015819927991942456,
     >  0.01597739164815558,
     >  0.0186283068226432,
     >  0.01597646679124881,
     >  0.0175694909521229,
     >  0.018325278551856783,
     >  0.009333750247599235,
     >  0.015077324141077597,
     >  0.014742688264559281,
     >  0.014628337002312931,
     >  0.01825121447048229,
     >  0.014576154497259297,
     >  0.014983169037060444,
     >  0.017427816480509742,
     >  0.01614850357517326,
     >  0.01523786026187711,
     >  0.020085798161947218,
     >  0.017464325481998416,
     >  0.011038500454808779,
     >  0.016151195946048497,
     >  0.019202429699250914,
     >  0.01583952767232326,
     >  0.017731422595169355,
     >  0.012216914175789573,
     >  0.019746928474867184,
     >  0.02090919845720336,
     >  0.015581614746626643,
     >  0.0168801004757302,
     >  0.01781539632193191,
     >  0.014624524644823972,
     >  0.01681024088933729,
     >  0.021933786020444504,
     >  0.015863471254502992,
     >  0.013303421722986278,
     >  0.017815462370066468,
     >  0.021420603939579102,
     >  0.023020095370882343,
     >  0.01948515170158825,
     >  0.022430367699889455,
     >  0.015417378283528528,
     >  0.02110590609950101,
     >  0.015948270678408094,
     >  0.02045435798468151,
     >  0.02013657050325178,
     >  0.015968626148256837,
     >  0.01597737922139264,
     >  0.018081928602444827,
     >  0.0158474107154917,
     >  0.019685771388903237,
     >  0.016116960005054304,
     >  0.0159633356693918,
     >  0.021380308727861855,
     >  0.015393985773032609,
     >  0.02004389251768653,
     >  0.01882938383816643,
     >  0.015260929161675691,
     >  0.015681028949560652,
     >  0.01535110613914116,
     >  0.01617327379447871,
     >  0.014923522590865088,
     >  0.01790750284382865,
     >  0.017848894085455175,
     >  0.015582287615343457,
     >  0.012246226899244485,
     >  0.017859978456408132,
     >  0.014856564208885994,
     >  0.020028117946321842,
     >  0.02021737759474557,
     >  0.015512195816817025,
     >  0.01782767305594357,
     >  0.017429473804186524,
     >  0.012999172131635732,
     >  0.0175916082743183,
     >  0.016146681922181202,
     >  0.017592138162250387,
     >  0.017709490364715544,
     >  0.015563517192083487,
     >  0.021455870517361437,
     >  0.017424243976076004,
     >  0.015977393169022816,
     >  0.019202445391392837,
     >  0.01862264949588148,
     >  0.013627903790051174,
     >  0.01785402964110412,
     >  0.013159084118914005
     >     /

      data bNs /
     >   0.014434178288728554,
     >   0.015296621431664291,
     >   0.014553612174587116,
     >   0.01471385435381519,
     >   0.013030680541774343,
     >   0.015740509327736872,
     >   0.0148944352383832,
     >   0.014868894886400509,
     >   0.015352918029181144,
     >   0.012573283834648612,
     >   0.014650113180923616,
     >   0.015565639410673475,
     >   0.013037399492995358,
     >   0.015489187523990575,
     >   0.013474207881884442,
     >   0.015220778162131791,
     >   0.01410491753383031,
     >   0.014593161499203075,
     >   0.0129121028246414,
     >   0.01483405853132991,
     >   0.01439031012971474,
     >   0.011993780548751843,
     >   0.012996666886038255,
     >   0.014014725067695771,
     >   0.015015855306125698,
     >   0.013490901549037487,
     >   0.016253793153007998,
     >   0.01562898716456286,
     >   0.0130871924912302,
     >   0.01583991155276586,
     >   0.013757132848923631,
     >   0.014653720117648041,
     >   0.014372176242566032,
     >   0.013917280575896758,
     >   0.013581271992092012,
     >   0.015126499552446547,
     >   0.016280339837122154,
     >   0.011721824655364914,
     >   0.013956742245538029,
     >   0.014181955521228471,
     >   0.015755178635016017,
     >   0.01355980915597655,
     >   0.014495947959615898,
     >   0.014044662259364286,
     >   0.014601136779860828,
     >   0.014358279281757225,
     >   0.016742064256423838,
     >   0.014390755973871459,
     >   0.015980183063361117,
     >   0.015020575577401218,
     >   0.017083955179459438,
     >   0.013122866771703469,
     >   0.014062700301586468,
     >   0.01262794022959237,
     >   0.009400969492673423,
     >   0.014573390457934944,
     >   0.010668892260888685,
     >   0.013564102509053523,
     >   0.013946391508881918,
     >   0.01498594678409917,
     >   0.013252572535934072,
     >   0.013958848609410338,
     >   0.015971968294982257,
     >   0.012261983127166032,
     >   0.01630237869654944,
     >   0.014425844022237566,
     >   0.01286618800977922,
     >   0.01516906999721688,
     >   0.014353223775944128,
     >   0.01304965704664668,
     >   0.013619713004166922,
     >   0.015150996293568109,
     >   0.012344116449233831,
     >   0.01434020572810224,
     >   0.01611787208803617,
     >   0.013063050309010706,
     >   0.014502835857967527,
     >   0.015758658396773206,
     >   0.012848137670415484,
     >   0.013361821242721427,
     >   0.015280381655405906,
     >   0.01751188678077182,
     >   0.013981521824946823,
     >   0.01507483006295294,
     >   0.015225571496496941,
     >   0.015681258010830472,
     >   0.014854877777138809,
     >   0.015916745657716738,
     >   0.014660743104776872,
     >   0.013979177464263276,
     >   0.010882605576580677,
     >   0.013192459259680149,
     >   0.015238528949039347,
     >   0.012637249977846259,
     >   0.014141636726544697,
     >   0.01443255829670984,
     >   0.012830120198109578,
     >   0.01579733412259189,
     >   0.013174690779960578,
     >   0.015490147444713354,
     >   0.01427012186000409,
     >   0.01291707156311517,
     >   0.0154841475407351,
     >   0.012500621929837629,
     >   0.013587570757806408,
     >   0.013270232894235762,
     >   0.015774339769579662,
     >   0.01408383275442577,
     >   0.015306868216486476,
     >   0.015079745431314302,
     >   0.013066814366234731,
     >   0.010716308616003596,
     >   0.014755933568651925,
     >   0.01172661338055863,
     >   0.01243116463561407,
     >   0.014206441868178129,
     >   0.015382661972868687,
     >   0.013576064617411223,
     >   0.017584012366349613,
     >   0.014249190179494862,
     >   0.012691421336972672,
     >   0.013120517675441984,
     >   0.013653041951250472,
     >   0.016101135293493458,
     >   0.01376556160244398,
     >   0.014352966418827863,
     >   0.01255544695540229,
     >   0.015101605098056248,
     >   0.010871735199605293,
     >   0.014048790807350156,
     >   0.016129806534927658,
     >   0.013611538255213771,
     >   0.01452520289565788,
     >   0.013146055271329123,
     >   0.014923851309006747,
     >   0.011769415921019856,
     >   0.014836585505687663,
     >   0.014141765892462237,
     >   0.014963772961532312,
     >   0.017328754014077345,
     >   0.01220510674266715,
     >   0.014863498984600687,
     >   0.012221454960923077,
     >   0.015413466801203068,
     >   0.014540457133152553,
     >   0.013403320221150485,
     >   0.013715024959490271,
     >   0.014071870153715095,
     >   0.014321837305615231,
     >   0.015196451173994548,
     >   0.014733673825122124,
     >   0.01526612537897129,
     >   0.013758572620382587,
     >   0.015407722875296785,
     >   0.013759076032327532,
     >   0.015387362343422438,
     >   0.01568741133337571,
     >   0.014035165529867484,
     >   0.013296800931647855,
     >   0.013801807288640888,
     >   0.012412022949627144,
     >   0.014162283047505206,
     >   0.015074356505789276,
     >   0.0134910092686559,
     >   0.01577844023813524,
     >   0.015331908229340259,
     >   0.014833813194964053,
     >   0.014592767812078175,
     >   0.013264880713166102,
     >   0.015907115119334417,
     >   0.014240312972830453,
     >   0.014780213961671528,
     >   0.012543450767225005,
     >   0.013431347596407768,
     >   0.015922918073799468,
     >   0.013595442730908007,
     >   0.01406080479016364,
     >   0.014067414026875294,
     >   0.015368321208485252,
     >   0.015049675813449966,
     >   0.014675615397603231,
     >   0.015452644341939813,
     >   0.012881672193898248,
     >   0.014485352074738538,
     >   0.015564597156574475,
     >   0.0156465445144212,
     >   0.015021824120245041,
     >   0.01531194028105526,
     >   0.01538101157930158,
     >   0.015178342744934903,
     >   0.012939645437208769,
     >   0.01618995946663538,
     >   0.015454198873883203,
     >   0.011378571527913764,
     >   0.010496775789644723,
     >   0.014562263547368427,
     >   0.0145513944026266,
     >   0.013351988459499135,
     >   0.016890263061731914,
     >   0.014499303157080472,
     >   0.015328598762916574
     >         /

      nloops = 2
      nll = 3
      prepdf = 1
      preff  = 1

      CALL InitPDFsetByNameM(1,"CT14nlo")
      CALL InitPDFM(1,0)
      CALL InitPDFsetByNameM(2,"EPPS16HE")
      CALL InitPDFM(2,0)
      CALL InitPDFsetByNameM(3,"EPPS16NE")
      CALL InitPDFM(3,0)
      CALL InitPDFsetByNameM(4,"EPPS16KR")
      CALL InitPDFM(4,0)
      CALL InitPDFsetByNameM(5,"EPPS16XE")
      CALL InitPDFM(5,0)
      CALL InitPDFsetByNameM(6,"LIKEnHE")
      CALL InitPDFM(6,0)
      CALL InitPDFsetByNameM(7,"LIKEnNE")
      CALL InitPDFM(7,0)
      CALL InitPDFsetByNameM(8,"LIKEnKR")
      CALL InitPDFM(8,0)
      CALL InitPDFsetByNameM(9,"LIKEnXE")
      CALL InitPDFM(9,0)
      CALL InitPDFsetByNameM(10,"EPPS16BE")
      CALL InitPDFM(10,0)
      CALL InitPDFsetByNameM(11,"EPPS16FE")
      CALL InitPDFM(11,0)
      CALL InitPDFsetByNameM(12,"EPPS16WW")
      CALL InitPDFM(12,0)
      CALL InitPDFsetByNameM(13,"EPPS16JCC")
      CALL InitPDFM(13,0)
      CALL InitPDFsetByNameM(14,"EPPS16JFE")
      CALL InitPDFM(14,0)
      CALL InitPDFsetByNameM(15,"EPPS16JPB")
      CALL InitPDFM(15,0)
      CALL InitPDFsetByNameM(16,"LIKEnCC")
      CALL InitPDFM(16,0)
      CALL InitPDFsetByNameM(17,"LIKEnFE")
      CALL InitPDFM(17,0)
      CALL InitPDFsetByNameM(18,"LIKEnPB")
      CALL InitPDFM(18,0)
      CALL InitPDFsetByNameM(20,"EPPS16AU")
      CALL InitPDFM(20,0)
      CALL InitPDFsetByNameM(21,"EPPS16PR")
      CALL InitPDFM(21,0)
      CALL InitPDFsetByNameM(22,"EPPS16CA")
      CALL InitPDFM(22,0)
      CALL InitPDFsetByNameM(23,"LIKEnVC")
      CALL InitPDFM(23,0)
      CALL InitPDFsetByNameM(24,"LIKEnAU")
      CALL InitPDFM(24,0)

      zh = 0.05
      nloops = 1
      Q = 2.

      IH = 1
      IC = 1

            
      H_AA = 1
      APFEL_init = dsqrt(1d0) - eps

C--------- Configure APFEL:
      call SetTimeLikeEvolution(.true.)
      call SetVFNS
      call SetPerturbativeOrder(1)
      call SetPDFEvolution("exactalpha")
C----- The next line tells APFEL to look for a function called "ExternalSetAPFEL"
      call SetPDFSet("external") !
      call SetMaxFlavourPDFs(5)
      call SetMaxFlavourAlpha(5)
      call InitializeAPFEL

      call EvolveAPFEL(APFEL_init ,Q)
      call  Nuclear_TMDFF(zh,Q,H_AA,fu,fub,fd,fdb,fs,fsb,fg,nloops)

      print *, "Apfel", fu

      call FF(zh,Q,fu,fub,fd,fdb,fs,fsb,fg,nloops)

      print *, "DSS14", fu

      call LIKEn(1,1,1,0,zh,Q**2.,1,fU,fUB,fD,fDB,fS,fSB,
     >            fC,fCB,fB,fBB,fG)

      print *, "LIKEn", fU/zh

      H_AA = 1

C--------- Configure APFEL:
      call SetTimeLikeEvolution(.true.)
      call SetVFNS
      call SetPerturbativeOrder(1)
      call SetPDFEvolution("exactalpha")
C----- The next line tells APFEL to look for a function called "ExternalSetAPFEL"
      call SetPDFSet("external") !
      call SetMaxFlavourPDFs(5)
      call SetMaxFlavourAlpha(5)
      call InitializeAPFEL

      call EvolveAPFEL(APFEL_init ,Q)
      call  Nuclear_TMDFF(zh,Q,H_AA,fu,fub,fd,fdb,fs,fsb,fg,nloops)

C      print *, fu

      END PROGRAM EIC
