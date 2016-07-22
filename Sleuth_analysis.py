#!/usr/bin/env python
import sys,os
import csv
import matplotlib.pyplot as plt
import math
import pylab as P
import numpy as np
import glob


#Rese@rch44_4me

def sort_table(table, col, backward):
    return sorted(table, key=operator.itemgetter(col), reverse=backward)

def columntypes(row):
    global col_names
    result=[]

    #row.pop()
    for column in range(len(row)):
        try:
            val=float(row[column])

        except:
            val=row[column]
            #print "column not converted", column
        
        #if val > 200000.0:
            #val=200000.0
        if val == "":
            val = 0.0
        if isinstance(val,str):
            if val == "NA":
                if (column == col_names["padj_AB"]) or (column == col_names["padj_AC"]) or (column == col_names["padj_BC"]): 
                    val = 1.0
                if (column == col_names["log2FoldChange_AB"]) or (column == col_names["log2FoldChange_AC"]) or (column == col_names["log2FoldChange_BC"]): 
                    val = 0.0
        #elif val > 200000.0:
            #val=200000.0
            
        if column == 1:
            result.append(row[column])
        else:
            result.append(val)
        
    pAB=result[col_names["padj_AB"]]
    pAC=result[col_names["padj_AC"]]
    pBC=result[col_names["padj_BC"]]
    fAB=result[col_names["log2FoldChange_AB"]]
    fAC=result[col_names["log2FoldChange_AC"]]
    fBC=result[col_names["log2FoldChange_BC"]]
    
    all_vals=result[2:11]

    root="DKAG00"
    exp=["A","B","C"]
    repnum=["1","2","3"]
    fudge=0.0
    for cond in exp:
        reps=[]
        originals=[]
        for rep in repnum:
            reps.append(result[col_names[root+rep+cond]]+fudge)
            originals.append(result[col_names[root+rep+cond]])
            
        avg=np.mean(reps)
        std=np.std(reps)
        med=np.median(reps)
        if avg > 1:
            remove=False
            rpos=0
            biggest=0
            
            for x in range(len(reps)):
                spot=reps[x]
                dist=math.fabs(spot-avg)
                dist2=math.fabs(spot-med)
                #cut=2.0*std
                grub=dist/(std+0.1)
                #if spot > 100000.0:
                    #print "Why not remove?", spot, "From ", all_vals, grub
                if grub > 1.4:
                    remove=True
                    if grub > biggest:
                        biggest=grub
                        rpos=x
                    
                    #foldc=spot/(med)
                    #if foldc > 10:
                        
                        #print int(foldc),result[1],"Experiment", cond, "Cutting ", spot, x, "From ", reps, std
            if remove:
                bavg=np.mean(all_vals)
                std=np.std(all_vals)
                med=np.median(all_vals)
                if med < 1:
                    med=1.0
                dist=math.fabs(originals[rpos]-bavg)
                grub=dist/(std+0.1)
                foldc=originals[rpos]/(med)
                
                if grub > 2.5: 
                    #print int(foldc),result[1],"Experiment", cond, "Cutting ", int(originals[rpos]), "From ", all_vals
                    if rpos > 0:
                        new=originals[:rpos]
                    else:
                        new=[]
                    rest=originals[rpos+1:]
                    if len(rest) > 0:
                        for thing in rest:
                            new.append(thing)
                    #print "new = ", new
                    avg=np.mean(new)
                    replace=str(rpos+1)
                    result[col_names[root+replace+cond]]=avg
        if avg < 1:
            ans=0
        else:
            ans=math.log(avg,2)
        result.append(ans)
    
    
    direc=0
    if pAB < 0.05:
        if fAB > 1.0:
            direc=1
        if fAB < -1.0:
            direc=-1
    result.append(direc)
    
    direc=0
    if pAC < 0.05:
        if fAC > 1.0:
            direc=1
        if fAC < -1.0:
            direc=-1
    result.append(direc)
    
    direc=0
    if pBC < 0.05:
        if fBC > 1.0:
            direc=1
        if fBC < -1.0:
            direc=-1
    result.append(direc)
    
    return result

def add2dict(row):
    global raw_gene_dict, out_head, repname
    
    result=[]
    for column in range(len(row)):
        try:
            val=round(float(row[column]),4)
            #val=round(val,4)

        except:
            val=row[column]
            #print "column not converted", column
        
        result.append(val)
    #print result

    ID=result[0]
    #print ID
    #ID=str(ID[:-2])
    #ID=int(ID[-6:])
    #print ID
    
    count=result[3]
    tpm=result[4]
    
    if ID in raw_gene_dict.keys():
        current=list(raw_gene_dict[ID])
        current.append(count)
        current.append(tpm)
        raw_gene_dict[ID]=list(current)
        
    else:
        raw_gene_dict[ID]=[count,tpm]
    
    return result

def add2diff(row):
    global diff_dict
    
    result=[]
    for column in range(len(row)):
        try:
            val=round(float(row[column]),4)
            #val=round(val,4)

        except:
            val=row[column]
            #print "column not converted", column
        
        result.append(val)
    #print result

    ID=result[1]
    #print ID
    #ID=str(ID[:-2])
    #ID=int(ID[-6:])
    #print ID
    
    qval=result[3]
    mean=result[6]
    b=result[4]
    
    if ID in diff_dict.keys():
        current=list(diff_dict[ID])
        current.append(mean)
        current.append(qval)
        current.append(b)
        diff_dict[ID]=list(current)
        
    else:
        diff_dict[ID]=[mean,qval,b]
    
    return result

def removespace(txtlist):
    
    fixed=[]
    for item in txtlist:
        fixed.append(item.replace(" ","_"))
    
    return fixed

def getcolumn(subgroup,column):
    global col_names
    #column=col_names[columname]
    col=[]
    for row in subgroup:
        col.append(row[column])
    return col

def matreader(stream):
    global col_names, Header
    "Read in raw data table"
    #global Header, TFs, col_names, outf
    #Read in the gene
    #Should I make a dictionary of columns and headings? Probably
    reader=csv.reader(stream,delimiter='\t')
    reader= list(reader)
    #fix spaces...
    Header=reader[0]
    Header=removespace(Header)
    Header.append("A_knockdown_avg")
    Header.append("B_control_avg")
    Header.append("C_treatment_avg")
    Header.append("dir_AB")
    Header.append("dir_AC")
    Header.append("dir_BC")
    for x in range(len(Header)):
        col_names[Header[x]]=x
    #Header.pop()

    reader=reader[1:]

    reader=map(columntypes,reader)

    return reader

def abundancereader(stream):
    #global col_names, Header
    "Read in raw data table"
    reader=csv.reader(stream,delimiter='\t')
    reader= list(reader)
    
    reader=reader[1:]

    reader=map(add2dict,reader)

    return 0

def diffreader(stream):
    #global col_names, Header
    "Read in raw data table"
    reader=csv.reader(stream,delimiter='\t')
    reader= list(reader)
    
    reader=reader[1:]

    reader=map(add2diff,reader)

    return 0

def glistreader(stream):
    global glist_dict
    "Read in raw data table"
    reader=csv.reader(stream,delimiter='\t')
    reader= list(reader)
    
    reader=reader[1:]
    
    for line in reader:
        tid=str(line[1])
        #tid=int(tid[-6:])
        gname=line[2]
        #print tid, gname
        """
        #Apparently there a bunch of duplicates in here...
        if tid in glist_dict.keys():
            print "WTF!", tid, gname
        else:
            glist_dict[tid]=gname
            """
        glist_dict[tid]=gname
    #reader=map(add2dict,reader)

    return 0

def filterconds(column, condition, value):
    global col_names
    #column=col_names[columname]
    def removerows(row):
        if condition == "<":
            return row[column]<value
        elif condition == ">":
            return row[column]>value
        elif condition == "=":
            return row[column]== value
        elif condition == "not":
            return row[column] != value
        elif condition == "contains":
            return value in row[column]
        elif condition == "True":
            return row[column]
        elif condition == "False":
            return not row[column]
        else:
            print "Bad condition"
            return
        
    return removerows

def bestfit(X,Y):
    if len(X) != len(Y):
        print "Lists aren't same length, moron!!!"
        return 0.0, 0.0, 0.0
    
    N=len(X)
    sumX=0.0
    sumY=0.0
    sumXsq=0.0
    sumYsq=0.0
    sumXY=0.0
    
    for pos in range(N):
        x=float(X[pos])
        y=float(Y[pos])
        sumX=sumX+x
        sumY=sumY+y
        sumXsq=sumXsq+(x*x)
        sumYsq=sumYsq+(y*y)
        sumXY=sumXY+(x*y)
    
    n=float(N)    
    slope=(n*sumXY-sumX*sumY)/(n*sumXsq-sumX**2)
    intercept=(sumY-slope*sumX)/n
    R=(n*sumXY-sumX*sumY)/(math.sqrt(n*sumXsq-sumX**2)*math.sqrt(n*sumYsq-sumY**2))
    Rsq=R*R
    Rsq=round(Rsq,3)
    
    return slope, intercept, Rsq

def plotstraightline(X,slope,intercept):
    linex=[min(X), 0.0, max(X)]
    liney=[(linex[0]*slope+intercept), intercept, (linex[2]*slope+intercept)]
    return linex, liney

def putoutall(out):
    def putout(row):
        "write a list of lists to a file"
        for n in range(len(row)):
            if n < (len(row)-1):
                out.write(str(row[n])+"\t")
            else:
                out.write(str(row[n]))
        out.write("\n")
    return putout

def putoutrow(out,row):
    "write a list of lists to a file"
    for n in range(len(row)):
        if n < (len(row)-1):
            out.write(str(row[n])+"\t")
        else:
            out.write(str(row[n]))
    out.write("\n")
    
def putoutpart(out,row,max):
    "write a list of lists to a file"
    for n in range(max):
        out.write(str(row[n])+"\t")

def make_corr_plot(data):
    global Header
    
    num_var=0
    
    line=data[1]
    data_cols=[]
    for x in range(len(Header)):
        if "DK" in Header[x]:
            data_cols.append(x)
            
    
    #data_cols=data_cols[:5]
    
    num_var=len(data_cols)
    
    num_plots=0
    
    x=num_var-1
    while x > 0:
        num_plots=num_plots+x
        x=x-1
        
    dim=int(math.sqrt(num_plots))
    print num_var,"variables = ",num_plots,"plots",dim
    
    #cons=filter(filterconds(1, "contains", "C"), data)
    #stresses=filter(filterconds(1, "contains", "S"), data)
    
    fig = plt.figure(figsize=(10,10))
    num_done=1
    
    for x in range(num_var):
        do_con=True
        #do_stress=True
        column=data_cols[x]
        conxdata=getcolumn(data,column)
        #stressxdata=getcolumn(stresses,column)
        if sum(conxdata)<1:
            do_con=False
        #if sum(stressxdata)<1:
            #do_stress=False
        for y in range(x+1,num_var):
            column=data_cols[y]
            conydata=getcolumn(data,column)
            #stressydata=getcolumn(stresses,column)
            if sum(conydata)<1:
                do_con=False
            #if sum(stressydata)<1:
                #do_stress=False
            if do_con:
                #ax = fig.add_subplot(dim+1,dim,num_done)
                ax = fig.add_subplot(1,1,1)
                num_done=num_done+1
            
            best=0
            if do_con:
                plt.plot(conxdata, conydata, 'bo', alpha=.9, label='Controls')
                slope, intercept, Rsq=bestfit(conxdata,conydata)
                linex, liney=plotstraightline(conxdata,slope,intercept)
                plt.plot(linex, liney, 'blue', label="R^2= "+str(Rsq))
                best=int(Rsq*100.0)
            """    
            if do_stress:   
                plt.plot(stressxdata, stressydata, 'ro', alpha=.9, label='Stressed')
                slope, intercept, Rsq=bestfit(stressxdata,stressydata)
                linex, liney=plotstraightline(stressxdata,slope,intercept)
                plt.plot(linex, liney, 'red', label="R^2= "+str(Rsq))
                if int(Rsq*100.0) > best:
                    best=int(Rsq*100.0)
            """    
            if do_con:
                plt.xlabel(Header[data_cols[x]])
                plt.ylabel(Header[data_cols[y]])
                plt.legend(loc=1)
                plt.savefig(str(best)+"_"+Header[data_cols[x]]+"_v_"+Header[data_cols[y]]+"_.pdf")
                plt.clf()

    #plt.savefig("Corr_plots.pdf")
    return 0

def make_diff_corr_plot(data):
    global col_names
    
    
    #controls=getcolumn(data, col_names["B_control_avg"])
    #knocks=getcolumn(data, col_names["A_knockdown_avg"])
    #treats=getcolumn(data, col_names["C_treatment_avg"])
    
    current=filter(filterconds(col_names["dir_AB"], "<", 0), data)
    
    upx=getcolumn(current, col_names["B_control_avg"])
    upy=getcolumn(current, col_names["A_knockdown_avg"])
    print "UPS", len(current), len(upx), len(upy)
    
    current=filter(filterconds(col_names["dir_AB"], "=", 0), data)
    x=getcolumn(current, col_names["B_control_avg"])
    y=getcolumn(current, col_names["A_knockdown_avg"])
    
    current=filter(filterconds(col_names["dir_AB"], ">", 0), data)

    dx=getcolumn(current, col_names["B_control_avg"])
    dy=getcolumn(current, col_names["A_knockdown_avg"])
    print "DOWNS", len(current), len(dx), len(dy)
    
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(1,1,1)
    
    plt.plot(x, y, 'ko', alpha=.2, label='Unresponsive')
    #slope, intercept, Rsq=bestfit(x,y)
    #linex, liney=plotstraightline(x,slope,intercept)
    #plt.plot(linex, liney, 'green', label="R^2= "+str(Rsq))
    #best=int(Rsq*100.0)
    plt.plot(upx, upy, 'ro', alpha=.9, label='Up regulated')
    plt.plot(dx, dy, 'bo', alpha=.9, label='Down regulated')
    
    plt.xlabel("B_control_avg")
    plt.ylabel("A_knockdown_avg")
    plt.legend(loc=2)
    plt.savefig("A_knockdown_avg_v_B_control_avg_.pdf")
    plt.clf()
    
    current=filter(filterconds(col_names["dir_BC"], ">", 0), data)
    
    upx=getcolumn(current, col_names["B_control_avg"])
    upy=getcolumn(current, col_names["C_treatment_avg"])
    print "UPS", len(current), len(upx), len(upy)
    
    current=filter(filterconds(col_names["dir_BC"], "=", 0), data)
    x=getcolumn(current, col_names["B_control_avg"])
    y=getcolumn(current, col_names["C_treatment_avg"])
    
    current=filter(filterconds(col_names["dir_BC"], "<", 0), data)

    dx=getcolumn(current, col_names["B_control_avg"])
    dy=getcolumn(current, col_names["C_treatment_avg"])
    print "DOWNS", len(current), len(dx), len(dy)
    
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(1,1,1)
    
    plt.plot(x, y, 'ko', alpha=.2, label='Unresponsive')
    #slope, intercept, Rsq=bestfit(x,y)
    #linex, liney=plotstraightline(x,slope,intercept)
    #plt.plot(linex, liney, 'green', label="R^2= "+str(Rsq))
    #best=int(Rsq*100.0)
    plt.plot(upx, upy, 'ro', alpha=.9, label='Up regulated')
    plt.plot(dx, dy, 'bo', alpha=.9, label='Down regulated')
    
    plt.xlabel("B_control_avg")
    plt.ylabel("C_treatment_avg")
    plt.legend(loc=2)
    plt.savefig("C_treatment_avg_v_B_control_avg_.pdf")
    plt.clf()
    
    current=filter(filterconds(col_names["dir_AC"], ">", 0), data)
    
    upx=getcolumn(current, col_names["A_knockdown_avg"])
    upy=getcolumn(current, col_names["C_treatment_avg"])
    print "UPS", len(current), len(upx), len(upy)
    
    current=filter(filterconds(col_names["dir_AB"], "=", 0), data)
    x=getcolumn(current, col_names["A_knockdown_avg"])
    y=getcolumn(current, col_names["C_treatment_avg"])
    
    current=filter(filterconds(col_names["dir_AB"], "<", 0), data)

    dx=getcolumn(current, col_names["A_knockdown_avg"])
    dy=getcolumn(current, col_names["C_treatment_avg"])
    print "DOWNS", len(current), len(dx), len(dy)
    
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(1,1,1)
    
    plt.plot(x, y, 'ko', alpha=.2, label='Unresponsive')
    slope, intercept, Rsq=bestfit(x,y)
    linex, liney=plotstraightline(x,slope,intercept)
    plt.plot(linex, liney, 'green', label="R^2= "+str(Rsq))
    best=int(Rsq*100.0)
    plt.plot(upx, upy, 'ro', alpha=.9, label='Up regulated')
    plt.plot(dx, dy, 'bo', alpha=.9, label='Down regulated')
    
    plt.xlabel("A_knockdown_avg")
    plt.ylabel("C_treatment_avg")
    plt.legend(loc=2)
    plt.savefig("A_knockdown_avg_v_C_treatment_avg_.pdf")
    plt.clf()
    
    
    return 0

def build_raw_dict(mat):
    global raw_gene_dict, col_names, Header
    
    exmp=0
    
    
    for row in mat:
        gene_name=row[1]
        if "-" not in gene_name:
            #break
        #else:
            entry=list(row)

            
            if gene_name in raw_gene_dict.keys():
                this=[]
                this.append(list(raw_gene_dict[gene_name]))
                this.append(list(entry))
                raw_gene_dict[gene_name]=this
                #if exmp < 10:
                    #print raw_gene_dict[gene_name]
                    #exmp=exmp+1
            else:
                raw_gene_dict[gene_name]=entry
    
    return 0

def makecombohist(columnname,matrix):
    global col_names
    
    do_norm=True
    raw=getcolumn(matrix,col_names[columnname])
    print "Found this many for Histogram", len(raw)
    try:
        normed=getcolumn(matrix,col_names[columnname+"_norm"])
    except:
        do_norm=False
    #print len(normed)," = ", len(raw)
    
    """
    fig = plt.figure(figsize=(10,6))
    ax1 = fig.add_subplot(111)
    
    n, bins, patches = P.hist(raw, bins=100, normed=False, facecolor='b', alpha=0.5)
    if do_norm:
        n, bins, patches = P.hist(normed, bins=100, normed=False, facecolor='r', alpha=0.5)
    
    P.xlabel('expression level')
    P.ylabel('Genes')
    ax1.set_ylim(0,5000)
    P.legend()
    
    P.savefig(columnname+"_combo_hist.pdf")
    P.clf()
    """
    
    fig = plt.figure(figsize=(10,6))
    #f,(ax,ax2) = plt.subplots(2,1,sharex=True)
    ax1 = fig.add_subplot(111)
    
    n, bins, patches = P.hist(raw, bins=100, normed=False, facecolor='b', alpha=0.5)
    #print "Returned", n, bins, patches
    top=n[0]
    rest=n[1]
    
    P.xlabel('expression level')
    P.ylabel('Genes')
    ax1.set_ylim(0,3000)
    #ax.set_ylim(.9*tops,1.1*top)
    #ax2.set_ylim(0,1.2*rest)
    P.legend()
    
    P.savefig(columnname+"_hist.pdf")
    P.clf()
    
    if do_norm:
        fig = plt.figure(figsize=(10,6))
        ax1 = fig.add_subplot(111)
        
        n, bins, patches = P.hist(normed, bins=100, normed=False, facecolor='r', alpha=0.5)
        
        P.xlabel('expression level')
        P.ylabel('Genes')
        ax1.set_ylim(0,3000)
        P.legend()
    
    P.savefig(columnname+"_normed_hist.pdf")
    P.clf()
    
    return 0

####################################################START HERE!!!!!!!###############################################
try:
    gf = open("Mouse_genes.txt",'rU')
    #rgf = open("full_rat_list.txt",'rU')

    outf = open("Mike_compiled_small.txt",'w')
    #conf = open("condensed.txt",'w')
    

except:
    print'Cannot open file\n'
    sys.exit(0)

##Define some variables
Header=[] #List
out_head=["Gene","Transcript"]

col_names={} #Dictionary
raw_gene_dict={}
glist_dict={}
diff_dict={}

glistreader(gf)

names = glob.glob('./Abundances/*.tsv')
print "Reading Abundances"
for fname in names:
    #chrom,sep,rest=mess.partition(":")
    first,sep,rest=fname.partition("\\")
    repname,sep,junk=rest.partition("_abundance")
    print repname
    out_head.append(repname+"_count")
    out_head.append(repname+"_tpm")
    print fname
    g=open(fname)
    abundancereader(g)
    g.close()
    
    print "RAW", len(raw_gene_dict.keys())
    
difnames = glob.glob('*_dif.txt')
print "Reading Diffs"
for fname in difnames:

    repname,sep,junk=fname.partition("_")
    print repname
    
    out_head.append(repname+"_mean")
    out_head.append(repname+"_qval")
    out_head.append(repname+"_beta")
    g=open(fname)
    diffreader(g)
    g.close()

#print "Glist", glist_dict.keys()

#print "RAW", len(raw_gene_dict.keys())

putoutrow(outf, out_head)
    
for tid in raw_gene_dict.keys():
    
    if tid in glist_dict.keys():
        gname=glist_dict[tid]
    else:
        gname="missing"
        #print tid
    line=[]
    line.append(gname)
    line.append(tid)
    vals=list(raw_gene_dict[tid])
    for thing in vals:
        line.append(thing)
    
    if tid in diff_dict.keys():
        extra=list(diff_dict[tid])
    else:
        extra=["Missing"]
        
    for thing in extra:
        line.append(thing)    
        
    putoutrow(outf, line)
"""
#Read in data
data=matreader(inf)
build_raw_dict(data)
#make_corr_plot(data)
make_diff_corr_plot(data)
print len(raw_gene_dict)

putoutrow(outf,Header)
putoutrow(conf,Header)

map(putoutall(outf),data)

#print data[0]
#current=filter(filterconds(1, "not", "-"), data)
#print "Left with", len(current)
#make_corr_plot(data)

singles=0
multiples=0

for col in Header:
    if ("DK" in col) and ("_norm" not in col):
        #makecombohist(col,data)
        print col
        
    
for gene in raw_gene_dict.keys():
    this=raw_gene_dict[gene]
    #print len(this)
    if len(this) > 1 and len(this) < 25:
        multiples=multiples+1
        #print gene,len(this)
    else:
        putoutrow(conf,this)
        singles=singles+1


print singles, multiples
"""
