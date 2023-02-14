# Welcome to SLAdesigner
# Samuel Beppler
# 23 March 2020
# sla.py v1.1
# import argparse
import random
import Bio.Seq
from Bio.Seq import Seq
from Bio.Seq import MutableSeq as mSeq
from Bio.SeqUtils import MeltingTemp as TM
from primer3 import bindings as primer3


def sla(guidestrand: str) -> list:
    #########################################################################
    #########################################################################
    # P = "123456"
    # NTS = list("ATGC")
    STEMLOOP = "GTCGTATCCAGTGCAGGGTCCGAGGTATTCGCACTGGATACGAC"
    # ^ TODO Kramer 2011, validate sequence
    UR = "CCAGTGCAGGGTCCGAGGTA"  # TODO UR-20 Kramer 2011
    FLAP = "AATAAATCATAA"
    # ALT_FLAP = "CGGCGGC"    # The alt flap is good for increasing GC content and avoidng G-tetrads
    FLAP_N = 6

    gs = Seq(guidestrand).back_transcribe()
    # p = Seq(P)
    sl = Seq(STEMLOOP)

    rando = False

    # if (args.flap): FLAP = args.flap
    # if (args.alt_flap):
    #     FLAP = ALT_FLAP
    #     FLAP_N = 7

    def RTPrimer():
        # if args.rt:
        #     rt_n = args.rt
        # else:
        #     rt_n = 6  # from Kramer (2011)
        rt_n = 6  # from Kramer (2011)
        rtp = sl + gs.reverse_complement()[:rt_n]
        rtp_b = sl + gs.reverse_complement()[:rt_n + 1]
        rtp_c = sl + gs.reverse_complement()[:rt_n + 2]
        return rtp, rt_n, rtp_b, rtp_c

    def ForwardPrimer():
        # global rando
        # fp_14 = ""
        # fp_16 = ""
        fp_n = 12  # initial Forward Primer-cDNA overlap length = 12
        fp = Seq("")
        ldr = Seq(FLAP)
        score = 0
        flap_n = FLAP_N
        fp = ldr[-flap_n:] + gs[:fp_n]
        fp = mSeq(str(fp))
        # if (not args.alt_flap) and (TM.Tm_NN(fp, dnac1=250, dnac2=250) < 59):
        if TM.Tm_NN(fp, dnac1=250, dnac2=250) < 59:
            # Tm too low, extend and substitute GCs into 5' flap
            flap_n += 1
            fp.insert(0, ldr[-flap_n])
            i = flap_n - 1
            while TM.Tm_NN(fp, dnac1=250, dnac2=250) < 59 and i >= 0:
                # if args.flap: break  # Don't change 5' flap is user defined
                rando = True
                nt = random.choice('GC')
                fp.pop(i)
                fp.insert(i, nt)
                i -= 1
        # Extend FP to cover more of cDNA/GS sequence
        while TM.Tm_NN(fp, dnac1=250, dnac2=250) < 59 and fp_n < 17:
            fp_n += 1
            fp = fp[:flap_n] + gs[:fp_n]
        if TM.Tm_NN(fp, dnac1=250, dnac2=250) < 59:
            print("\t *****")
            print("WARNING! Low forward primer Tm.")
            print("Relax flap constraint or consider adjusting PCR conditions.")
            print("To view Tm use option -v")
            print("\t *****")
        elif TM.Tm_NN(fp, dnac1=250, dnac2=250) > 61:  # Tm too high, trim flap
            # This will never occur with ideal AT rich flap, even if all GC in GS
            while (TM.Tm_NN(fp, dnac1=250, dnac2=250) > 61) and (len(fp) > fp_n + 3):
                # Flap specified by user, just trim
                fp.pop(0)
                flap_n -= 1
            if TM.Tm_NN(fp, dnac1=250, dnac2=250) > 61:
                print("\t *****")
                print("WARNING! High forward primer Tm. Use default/AT rich flap.")
                print("Default 5' Flap: AATAAATCATAA")
                print("\t *****")
        # Extend leader sequence if within bound
        if ((flap_n < 7)
                and (TM.Tm_NN(Seq(ldr[-flap_n] + str(fp)), dnac1=250, dnac2=250) < 61)):
            fp.insert(0, ldr[-flap_n])
            flap_n += 1
        # Check for primer dimers
        homo_dg = primer3.calcHomodimer(str(fp)).dg  # delta-G in cal/mol
        het_dg = primer3.calcHeterodimer(str(fp), UR).dg
        # -9 kcal/mol is threshold delta-G (IDT)
        if (homo_dg <= -9000) or (het_dg <= -9000):
            # if (args.verbose): print("Dimer found. Shuffling flap.")
            k = 100
            while k > 0:
                # Shuffle 5' flap
                flap = list(str(fp[:-fp_n]))
                random.shuffle(flap)
                shuffled = mSeq("".join(flap)) + fp[-fp_n:]
                shuffled_homo_dg = primer3.calcHomodimer(str(shuffled)).dg
                shuffled_het_dg = primer3.calcHeterodimer(str(shuffled), UR).dg
                if min(shuffled_homo_dg, shuffled_het_dg) > min(homo_dg, het_dg):
                    fp = shuffled
                    homo_dg = shuffled_homo_dg
                    het_dg = shuffled_het_dg
                k -= 1
        if (fp_n == 12) or (fp_n == 13):  # Select experimental test components
            fp_2 = fp[:flap_n] + gs[:(fp_n + 2)]
            fp_3 = fp[:flap_n] + gs[:(fp_n + 4)]
        elif fp_n < 17:
            fp_2 = fp[:flap_n] + gs[:(fp_n - 2)]
            fp_3 = fp[:flap_n] + gs[:(fp_n + 2)]
        else:
            fp_2 = fp[:flap_n] + gs[:(fp_n - 3)]
            fp_3 = fp[:flap_n] + gs[:(fp_n - 1)]

        return fp.toseq(), fp_n, fp_2, fp_3, homo_dg, het_dg

    def Probe(fp, fp_n):
        # global gs
        # if args.probe:
        #     p_n = args.probe
        # else:
        #     p_n = len(gs[fp_n:])
        p_n = len(gs[fp_n:])
        i = 12 - p_n  # minimum length of stemloop to include
        p = mSeq(str(sl[-i:] + gs[-p_n:].reverse_complement()))
        # Ideally we would extend into stemloop until Tm is about 70 C
        # MGB data on exact Tm is limited, we estimate MGB by +20 C
        MGB_corr = 20
        # Limiting probe-stemloop orerlap length to 17 will avoid reverse primer (UR)
        while (i < 17) and (not TM.Tm_NN(p, dnac1=250, dnac2=250) + MGB_corr > 70):
            i += 1
            p.insert(0, sl[-i])
        if (len(gs) - fp_n) > 7:
            p_b = str(p)[:-2]
            p_b = sl[-(i + 2)] + sl[-(i + 1)] + p_b
        else:
            p_b = ""
        return p.toseq(), p_b

    rt, rt_n, rt_b, rt_c = RTPrimer()
    fp, fp_n, fp_b, fp_c, homo_dg, het_dg = ForwardPrimer()
    # Try k times if randomized flap
    if rando:
        # if args.k:
        #     k = args.k
        # else:
        #     k = 100
        k = 100
        for i in range(k):
            fp_k, fp_n_k, fp_b_k, fp_c_k, homo_dg_k, het_dg_k = ForwardPrimer()
            if min(homo_dg_k, het_dg_k) > min(homo_dg, het_dg):
                fp, fp_n, fp_b, fp_c, homo_dg, het_dg = fp_k, fp_n_k, fp_b_k, fp_c_k, homo_dg_k, het_dg_k
    p, p_b = Probe(fp, fp_n)
    FLAP_N = len(fp) - fp_n

    output = list()

    if (homo_dg <= -9000) or (het_dg <= -9000):
        temp = 'WARNING: Primer dimer found' + '\n'
        temp = temp + f'Homodimer delta-G (cal/mol): {homo_dg}\n'
        temp = temp + f'HHetero delta-G (cal/mol): {het_dg}'
        output.append(temp)

    output.append(f'RT-{rt_n}\t{rt}\t{TM.Tm_NN(rt):.2f}')
    output.append(f'RT-{rt_n + 1}\t{rt_b}\t{TM.Tm_NN(rt_b):.2f}')
    output.append(f'RT-{rt_n + 2}\t{rt_c}\t{TM.Tm_NN(rt_c):.2f}')
    output.append(f'F-{fp_n}\t{fp}\t{TM.Tm_NN(fp)}')
    output.append(f'F-{len(fp_b) - FLAP_N}\t{fp_b}\t{TM.Tm_NN(fp_b):.2f}')
    output.append(f'F-{len(fp_c) - FLAP_N}\t{fp_c}\t{TM.Tm_NN(fp_c):.2f}')
    try:
        output.append(f'P-{len(gs) - fp_n}\t{p}\t{TM.Tm_NN(p):.2f}')
    except IndexError as err:
        print(f'Unable to identify P-{len(gs) - fp_n}')
        output.append(f'P-{len(gs) - fp_n}\t{p}\t')
    try:
        output.append(f'P-{len(gs) - fp_n - 2}\t{p_b}\t{TM.Tm_NN(p_b):.2f}')
    except IndexError as err:
        print(f'Unable to identify P-{len(gs) - fp_n - 2}')
        output.append(f'P-{len(gs) - fp_n - 2}\t{p_b}\t')

    #########################################################################
    #########################################################################

    return output


def get_sla_output_numbers(sla_output: list) -> list:
    output = list()
    primers_out = sla_output[(len(sla_output) - 8):len(sla_output)]
    for i in range(len(primers_out)):
        temp = primers_out[i]
        if temp.endswith('\t'):
            output.append(0)
        else:
            start = temp.index('-') + len('-')
            end = temp.index('\t')
            num = int(temp[start:end])
            output.append(num)
    return output


def cross_check_sla(guidestrand: str, times=50) -> list:
    maxProbes = 0
    maxProbeLength = 0
    best_sla = list()
    for i in range(times):
        sla_out = sla(guidestrand)
        out_numbers = get_sla_output_numbers(sla_out)
        rt1, rt2, rt3, f1, f2, f3, p1, p2 = out_numbers
        num_probes = len(out_numbers) - out_numbers.count(0)
        probe_length = max(p1, p2)
        if num_probes == maxProbes:
            if probe_length > maxProbeLength:
                maxProbeLength = probe_length
                best_sla = sla_out
            # elif probe_length == maxProbeLength:
            #     pass
        elif num_probes > maxProbes:
            maxProbes = num_probes
            maxProbeLength = probe_length
            best_sla = sla_out
    return best_sla

