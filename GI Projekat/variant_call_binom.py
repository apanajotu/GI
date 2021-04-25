import sys
import numpy
import datetime
import re
import collections
from collections import OrderedDict
import math
import time

##
## STRING DEFINITIONS - constant strings of files that were used and additional info that couldn't be pulled from the files
##

ref_assembly_str = ',assembly=b37'
ref_md5_str = ',md5=bb77e60e9a492fd0172e2b11e6c16afd'
ref_species_str = ',species=\"Homo sapiens\"'
ref_taxonomy_str = ',taxonomy=x' #couldn't find what this acctually means


##
## VCF CONTIG INFO FUNC - Pulls out information on all contigs from reference .fasta.fai file
##

def write_contig_info(file):
    contig_info =''
    for line in file:
        field = line.split('\t')
        ID_str = field[0]
        len_str = field[1]
        contig_info += '##contig=<ID=' + ID_str + ",length=" + len_str +'>\n' # + ref_assembly_str + ref_md5_str + ref_species_str + ref_taxonomy_str 
    return contig_info
        
    
    
##
## VCF HEADER FUNCTION - Writes out header of the .vcf file
##

def write_vcf_header(file, fai_file_str):
    date_time = datetime.datetime.now()
    
    file_format_str = '##fileformat=VCF4.2\n'
    date_str = '##fileDate=' + date_time.strftime("%Y%m%d") + '\n'
    source_str = '##source=variant_call_binomV0.1\n'
    ref_file_str = '##reference=file://' + fai_file_str + '\n'
    contig_info_str = ''
    if(fai_file_str != ''):
        fai = open (fai_file_str, 'r')
        contig_info_str = write_contig_info(fai)
    else:
        print('Warning: Fasta Index file not provided. Contig info will not be available in .vcf file.')
    alt_str = '##ALT=<ID=*,Description="Represents allele(s) other than observed.">\n'
    indel_str = '##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Indicates that the variant is an INDEL.">\n'
    dp_str = '##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw read depth">\n'
    gt_str = '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    vaf_str = '##FORMAT=<ID=VAF,Number=1,Type=Float,Description="Variant frequency in sample">\n'
    
    table_header_str = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHCC1143BL\n'
    
    file.write(file_format_str + date_str + source_str + ref_file_str + contig_info_str + alt_str + indel_str + dp_str + gt_str + vaf_str + table_header_str)
    
    return file



##
## CALCULATE VARIANTS CALLS
##

class Variant:
    def __init__(self, symbol, num_occ, num_reads, ref_symbol):
        self.var_type, self.var_symbol, self.var_len = self.variant_def(symbol, ref_symbol)
        
        self.num_occ = num_occ
        self.num_reads = num_reads
        self.occ_percent = round(float(num_occ/num_reads), 1)
        self.VAF = str(self.occ_percent)
        self.DP = str(self.num_reads)
        
    def variant_def(self, symbol, ref_symbol):
        if ( len(symbol) == 1):
                var_symbol = symbol                
                if(symbol == '.'):
                    var_type = 'REF'
                else:
                    var_type = 'SNV'
                var_len = 1
        else:
            if(symbol[0]== '+'):
                var_type = 'INS'
            else:
                var_type = 'DEL'
            var_symbol = ref_symbol + ''.join(filter(str.isalpha, symbol))
            var_len = int(''.join(filter(str.isdigit, symbol)))
            
    
        return var_type, var_symbol, var_len
    

def binomial(n, k, p):
        return (math.factorial(n)/(math.factorial(k)*math.factorial(n-k)))*(p**k)*(1-p)**(n-k)
    

def calculate_binomial(var1, var2):
        res_vars = {}
        
        res1 = binomial(float(var1.num_reads), float(var1.num_occ), float(p))
        res2 = binomial(float(var2.num_reads), float(var2.num_occ), float(p))
        res3 = binomial(float(var1.num_reads), float(var1.num_reads), float(p))
        
        if(res1 > res2):
            if(res1 > res3):
                return var1, None
            else:
                return var1, var2
        else: #res1 < res2
            if(res2 > res3):
                return var2, None
            else:
                return var1, var2
    
def prepare_read_str(read_str): #raname to calculate_variants
        read_str = re.sub('\^.', '', read_str) #remove caret and next symbol
        read_str = re.sub('\$', '', read_str ) #remove dollar
        read_str = re.sub('\,', '.', read_str) #substitute comma with dot
        read_str = read_str.upper() #switch all letters to uppercase
        return read_str
    
##
## STORE SINGLE PILEUP CALLS
##


class Pileup_line:
    def __init__(self, new_line):
        self.new_line = new_line
        field = new_line.split("\t")
        self.seq_str = field[0] # sequence name string
        self.pos_str = field[1] # position string
        self.ref_str = field[2] # reference base string
        self.cnt_str = field[3] # read count string
        self.res_str = field[4] # result string
        self.qual_str = field[5] # base quality string
        
        self.pos = int(self.pos_str)
        self.cnt = int(self.cnt_str)
        
        self.read_vars = collections.OrderedDict()
        self.read_vars = {
            '.':0,
            'A':0,
            'T':0,
            'C':0,
            'G':0
        }
        
        self.var1 = None
        self.var2 = None
        
    def process_read_line(self):
        self.find_variants(self.res_str)
        return self.call_variants_binomial(self.var1, self.var2)
    
    
    def find_variants (self, read_str):
#         var1 = None
#         var2 = None
        has_ins = read_str.find('+')
        has_del = read_str.find('-')
        
        read_str = prepare_read_str(read_str)
        
        if ((has_ins == -1) and (has_del == -1)):            
            sym_cnt = collections.Counter(read_str.strip('"')).most_common()
#             print(sym_cnt)
            if (len(sym_cnt) == 1):
                self.var1 = Variant(sym_cnt[0][0], sym_cnt[0][1], self.cnt, self.ref_str)
                self.var2 = None
            else:
                self.var1 = Variant(sym_cnt[0][0], sym_cnt[0][1], self.cnt, self.ref_str)
                self.var2 = Variant(sym_cnt[1][0], sym_cnt[1][1], self.cnt, self.ref_str)
        
        elif (has_ins != -1 or has_del != -1):
#             print("INDEL")
            self.prepare_read_vars(read_str)
            if(self.read_vars_list[1][1] == 0):
                self.var1 = Variant(self.read_vars_list[0][0], self.read_vars_list[0][1], self.cnt, self.ref_str)
                self.var2 = None
            else:
                self.var1 = Variant(self.read_vars_list[0][0], self.read_vars_list[0][1], self.cnt, self.ref_str)
                self.var2 = Variant(self.read_vars_list[1][0], self.read_vars_list[1][1], self.cnt, self.ref_str)
                
#             print(self.var1.var_symbol)
#             print(self.var2.var_symbol)
                
                
        
    def prepare_read_vars(self, read_str):    
        skip_index = 0
        
        for i in range(0, len(read_str)):
            if(i<skip_index) : continue #skip indel symbols
             
            if(read_str[i] == '+' or read_str[i] == '-'):
                num_len = 0
                for j in range (i+1, len(read_str)):
                    if(read_str[j].isnumeric()): num_len +=1
                    else: break
                indel_len = int(read_str[i+1 : i+1+num_len])
#                 if(read_str[i-1] == '.'): 
                new_indel = read_str[(i):(i+num_len+indel_len+1)]
#                 else:
#                     new_ins = read_str[i-1] + read_str[(i+2):(i+2+ins_len)]
                if new_indel in self.read_vars : self.read_vars.update({new_indel:self.read_vars[new_indel]+1})
                else : self.read_vars.update({new_indel:1})
                self.read_vars.update({'.':self.read_vars['.']-1}) #need to substitute number of matches by one
                
                skip_index = i + num_len + indel_len +1
            
            if(read_str[i] == '.'): self.read_vars.update({'.' :self.read_vars['.']+1})
            if(read_str[i] == 'A'): self.read_vars.update({'A' :self.read_vars['A']+1})
            if(read_str[i] == 'T'): self.read_vars.update({'T' :self.read_vars['T']+1})
            if(read_str[i] == 'C'): self.read_vars.update({'C' :self.read_vars['C']+1})
            if(read_str[i] == 'G'): self.read_vars.update({'G' :self.read_vars['G']+1})
        
        self.read_vars = OrderedDict(sorted(self.read_vars.items(), key=lambda item: item[1], reverse=True))
        self.read_vars_list = list(self.read_vars.items())
        
    
    
    def call_variants_binomial(self, var1, var2):
        result_line = ''
        genotype_str = ''
        
        ref_field = 'ERROR'
        alt_field = 'ERROR'
        gt_field = 'ERROR'
        info_field = 'DP=' + var1.DP
        
        if var2 is None:   #single variant call      
#             print('single variant: ' + var1.var_symbol)
            ref_field = self.ref_str
            alt_field = var1.var_symbol
            if(var1.var_symbol == '.'):
                gt_field = '0/0:'+var1.VAF
            else:
                gt_field = '1/1:'+var1.VAF
                
                
                
        else:   #binomial calculation of two most probable calls     
            var1, var2 = calculate_binomial(var1, var2)
            
            if(var2 is None): #one call is the most probable one
#                 print('binomial, single variant: ' + var1.var_symbol)
                if(var1.var_type == 'INS' or var1.var_type == 'DEL'):
                    info_field = 'INDEL;' + info_field
                if(var1.var_type == 'DEL'): 
                    ref_field = var1.var_symbol
                    alt_field = self.ref_str
                    gt_field = '1/1:'+var1.VAF
                else:
                    ref_field = self.ref_str
                    alt_field = var1.var_symbol
                    if(var1.var_symbol == '.'):
                        gt_field = '0/0:'+var1.VAF
                    else:
                        gt_field = '1/1:'+var1.VAF
                        
            else: #two calls are most probable
#                 print('binomial, double variant: '+var1.var_symbol + '/' + var2.var_symbol)
                if(var1.var_type == 'INS' or var1.var_type == 'DEL' or var2.var_type == 'INS' or var2.var_type == 'DEL'):
                    info_field = 'INDEL;' + info_field
                    
                if(var1.var_type == 'DEL' and var2.var_type == 'DEL'):
                    ref_field = var1.var_symbol
                    alt_field = self.ref_str
                    gt_field = '1/1:'+var1.VAF
                    result_line += self.seq_str + '\t' + self.pos_str + '\t' + '.' + '\t' + ref_field + '\t' 
                    result_line += alt_field + '\t' + '.' + '\t' + 'PASS' + '\t' + info_field + '\t' + 'GT:VAF' + '\t' + gt_field +  '\n'
                    
                    ref_field = var2.var_symbol
                    alt_field = self.ref_str
                    gt_field = '1/1:'+var2.VAF
                    
                elif(var1.var_type == 'DEL' and var2.var_type != 'DEL'):
                    ref_field = var1.var_symbol
                    alt_field = self.ref_str
                    gt_field = '0/0:'+var1.VAF
                    result_line += self.seq_str + '\t' + self.pos_str + '\t' + '.' + '\t' + ref_field + '\t' 
                    result_line += alt_field + '\t' + '.' + '\t' + 'PASS' + '\t' + info_field + '\t' + 'GT:VAF' + '\t' + gt_field +  '\n'
                    
                    info_field = 'DP=' + var1.DP
                    ref_field = self.ref_str
                    alt_field = var2.var_symbol
                    if(var2.var_symbol == '.'):
                        gt_field = '0/0:' + var2.VAF
                    else:
                        gt_field = '0/1:' + var2.VAF
                
                elif(var1.var_type != 'DEL' and var2.var_type == 'DEL'):
                    ref_field = var2.var_symbol
                    alt_field = self.ref_str
                    gt_field = '0/0:'+var2.VAF
                    result_line += self.seq_str + '\t' + self.pos_str + '\t' + '.' + '\t' + ref_field + '\t' 
                    result_line += alt_field + '\t' + '.' + '\t' + 'PASS' + '\t' + info_field + '\t' + 'GT:VAF' + '\t' + gt_field +  '\n'
                    
                    info_field = 'DP=' + var1.DP
                    ref_field = self.ref_str
                    alt_field = var1.var_symbol
                    if(var1.var_symbol == '.'):
                        gt_field = '0/0:'+var1.VAF
                    else:
                        gt_field = '0/1:'+var1.VAF
                        
                else: # (var1.var_type != 'DEL' && var2.var_type != 'DEL')
                    ref_field = self.ref_str
                    if(var1.var_symbol == '.'):
                        alt_field = var2.var_symbol
                        gt_field = '0/1:' + var1.VAF + ',' + var2.VAF
                    elif (var2.var_symbol == '.'):
                        alt_field = var1.var_symbol
                        gt_field = '0/1:' + var2.VAF + ',' + var1.VAF
                    else:
                        alt_field = var1.var_symbol + ',' + var2.var_symbol
                        gt_field = '1/2:'+var1.VAF + ',' + var2.VAF
                            
                
        result_line += self.seq_str + '\t' + self.pos_str + '\t' + '.' + '\t' + ref_field + '\t' 
        result_line += alt_field + '\t' + '.' + '\t' + 'PASS' + '\t' + info_field + '\t' + 'GT:VAF' + '\t' + gt_field +  '\n'
        
        return result_line
    


##
## MAIN Function - goes through the pileup file, writes to .vcf file
## 
    
def main_func(pileup_file_str, output_file_str, reference_fai_file_str):
    with open (pileup_file_str, 'r') as pileup_file, open (output_file_str, 'w') as output_file:

        output_file = write_vcf_header(output_file, reference_fai_file_str)   
        line_num = 0;
        for line in iter(pileup_file.readline, ''):
            line_num+=1 
            #if (line_num > 3): break
            base_pileup = Pileup_line(line)
            if(base_pileup.cnt > 0):
                res_str = str(base_pileup.res_str)
                processed_str = base_pileup.process_read_line()
                output_file.write(processed_str)       

        print('Processing done.' + output_file_str + ' file created.')
        
        
        
##
## MAIN CODE - opens files, calls main task
##
if __name__ == "__main__" :
    if(len(sys.argv) < 2):
        print('Error: Please provide pileup file path')
        raise NameError('Missing_arguments')
    else:
        pileup_file_str = sys.argv[1]
        output_file_base_str = 'binom_variant'
        fai = ''
        p_array = [0.85]
        o_def = 0
        p_def = 0
        f_def = 0
        for i in range(2, len(sys.argv)):
            if(sys.argv[i-1] == '-o' or sys.argv[i-1] == '-f'): continue
            if(sys.argv[i] == '-o'):
                if(o_def == 1):
                        print('Error: -o switch already defined')
                        raise NameError('Incorrect_arguments')
                output_file_base_str = sys.argv[i+1]
                o_def = 1
            if(sys.argv[i] == '-f'):
                if(f_def == 1):
                        print('Error: -f switch already defined')
                        raise NameError('Incorrect_arguments')
                fai = sys.argv[i+1]
                f_def = 1
            if(sys.argv[i] == '-p'):
                if(p_def != 0):
                    print('Error: -p switch already defined')
                    raise NameError('Incorrect_arguments')
                p_def = i+1
                
        if(p_def != 0): 
            p_array = []    
            for i in range(p_def, len(sys.argv)):    
                if(sys.argv[i] == '-o' or sys.argv[i] == '-f'): break
                p_array.append(float(sys.argv[i]))
            if(len(p_array) == 0):
                print('Error: No percentages defined')
                raise NameError('Incorrect_arguments')
            
        start = time.time()
        for i in p_array:
            i_str = str(i)
            p = i
            if (output_file_base_str.rfind(".vcf", 0, len(output_file_base_str)) == -1):
                output_file_str = output_file_base_str + '_p' + str(i) + '_called.vcf'
            else:
                output_file_str = output_file_base_str[0:output_file_base_str.rfind(".vcf", 0, len(output_file_base_str))] + '_p' + str(i) + '_called.vcf'
            main_func(pileup_file_str, output_file_str, fai)
        end = time.time() 
        print('Processing took ' + str(end - start) + ' seconds')