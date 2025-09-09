import SelCarbon as sc
import del_overthreebridge as del3
import HSC as hsc
import pandas as pd

need_C_num = int(input("Please enter the required carbon number:"))
need_hetero = input("Whether heteroatom is required:y/n\n")

if need_hetero == 'n':
    df_smi = pd.DataFrame()
    df_bde = pd.DataFrame()
    for C_num in range(6,need_C_num+1):
        pd.set_option('display.max_columns', None)
        pd.set_option('display.width', 5000)
        Start = sc.SiftCarbonRings(C_num)
        CH_MF = Start.find_CH()
        carbon_rings = Start.choosen_CH(CH_MF)
        alfabetsmi = Start.smitoalfbetsmi(carbon_rings)
        pred_BDE = Start.carbon_bde(alfabetsmi)
        final = Start.none_rings_CH(pred_BDE)
        df_smi = pd.concat([df_smi, carbon_rings], axis=0)
        df_bde = pd.concat([df_bde, final], axis=0)
    delbridge = del3.del_overthreebridge(df_smi)
    hydsto = hsc.cal_Hweight(delbridge)
    hystoover6 = hsc.HSCover6(hydsto)
    screen_results = hsc.underMCH(hystoover6,df_bde)

elif need_hetero == 'y':
    heteroatom = input("Please enter the desired heteroatom:")
    df_smi = pd.DataFrame()
    df_bde = pd.DataFrame()
    for C_num in range(5,need_C_num+1):
        Start = sc.SiftCarbonHeteroatom(C_num,heteroatom)
        CH_MF = Start.find_C_hetero()
        carbon_rings = Start.choosen_heteroatom_ring(CH_MF)
        alfabetsmi = Start.smitoalfbetsmi(carbon_rings)
        pred_BDE = Start.carbon_bde(alfabetsmi)
        final = Start.none_rings_CH(pred_BDE)
        df_smi = pd.concat([df_smi, carbon_rings], axis=0)
        df_bde = pd.concat([df_bde, final], axis=0)
    delbridge = del3.del_overthreebridge(df_smi,heteroatom)
    hydsto = hsc.cal_Hweight(delbridge,heteroatom)
    hystoover6 = hsc.HSCover6(hydsto)
    screen_results = hsc.underMCH(hystoover6,df_bde)
