import time
import sys
import psutil

opt1 = []

class EfficientSequenceAlignment:
    def alignSequence(self,str_1,str_2):
        m = len(str_1)
        n = len(str_2)

        opt=[]

        opt = [[0] * (n + 1) for i in range(m + 1)]

        alpha = {'AA': 0, 'CC': 0, 'GG': 0, 'TT': 0, 
                'AC': 110, 'CA': 110, 'AG': 48, 'GA': 48, 
                'AT': 94, 'TA': 94, 'CG': 118, 'GC': 118,
                'CT': 48, 'TC': 48, 'GT': 110, 'TG': 110}
        delta = 30

        for i in range(1, m + 1):
            opt[i][0] = i*delta
        for j in range(0, n + 1):
            opt[0][j] = j*delta

        for i in range(1, m + 1):
            for j in range(1, n + 1):
                string = str_1[i - 1] + str_2[j - 1]
                opt[i][j] = min(opt[i - 1][j - 1] + alpha[string],
                               opt[i][j-1] + delta, opt[i-1][j] + delta)

        i, j = m, n
        x = ""
        y = ""

        while i and j:
            string = str_1[i - 1] + str_2[j - 1]
            if opt[i][j] == opt[i-1][j-1] + alpha[string]:
                x = str_1[i - 1] + x
                y = str_2[j - 1] + y
                i -= 1
                j -= 1
            elif opt[i][j] == opt[i - 1][j] + delta:
                x = str_1[i - 1] + x
                y = "_" + y
                i -= 1
            elif opt[i][j] == opt[i][j - 1] + delta:
                x = "_" + x
                y = str_2[j - 1] + y
                j -= 1

        while i:
            x = str_1[i - 1] + x
            y = "_" + y
            i -= 1

        while j:
            x = "_" + x
            y = str_2[j - 1] + y
            j -= 1

        return x, y, opt[m][n]
    
    def DC_sol(self, str_1, str_2):
        m = len(str_1)
        n = len(str_2)
        if m < 2 or n < 2:
            return self.alignSequence(str_1, str_2)
        else:
            left_part = self.align_space_efficient(str_1[:m // 2], str_2, 0)
            right_part = self.align_space_efficient(str_1[m // 2:], str_2, 1)
            
            new_list = [left_part[j] + right_part[n - j] for j in range(n + 1)]
            
            min_index = new_list.index(min(new_list))
            
            left_call = self.DC_sol(str_1[:len(str_1) // 2], str_2[:min_index])
            right_call = self.DC_sol(str_1[len(str_1) // 2:], str_2[min_index:])
            
            l = [left_call[r] + right_call[r] for r in range(3)]
            
        return [left_call[r] + right_call[r] for r in range(3)]


    def align_space_efficient(self, X, Y, flag):
        m = len(X)
        n = len(Y)
        I = 0

        alpha = {'AA': 0, 'CC': 0, 'GG': 0, 'TT': 0, 
                 'AC': 110, 'CA': 110, 'AG': 48, 'GA': 48, 
                 'AT': 94, 'TA': 94, 'CG': 118, 'GC': 118, 
                 'CT': 48, 'TC': 48, 'GT': 110, 'TG': 110}
        delta = 30
            
        for j in range(n + 1):
            opt1[0][j] = j*delta       
        
        if flag == 0:
            for i in range(1, m + 1):
                I=i%2
                opt1[I][0] = i*delta
                Iminus=(I-1)%2
                
                for j in range(1, n + 1):
                    opt1[I][j] = min(opt1[Iminus][j - 1] + alpha[X[i - 1] + Y[j - 1]],
                                   opt1[I][j - 1] + delta,
                                   opt1[Iminus][j] + delta)
                    
        elif flag == 1:
            for i in range(1, m + 1):
                
                I=i%2
                opt1[I][0] = i*delta
                Iminus=(I-1)%2
                
                for j in range(1, n + 1):
                    opt1[I][j] = min(opt1[Iminus][j - 1] + alpha[X[m - i] + Y[n - j]],
                                   opt1[I][j - 1] + delta,
                                   opt1[Iminus][j] + delta)

                
        return opt1[m%2]

    def input_generate(self, inputFileLocation):
        li = []

        with open(inputFileLocation) as f:
            lines = [x.rstrip("\n") for x in f.readlines()]
            for j in range(0, len(lines)):
                if(lines[j].isdigit()==False):
                    strind=j 
                else:
                    k = int(lines[j])
                    strr = lines[strind][0:k+1] + lines[strind] +  lines[strind][k+1:len(lines[strind][:])]        
                    lines[strind] = strr
            li.append(lines)

        li = li[0]
        f.close()
        return li


def process_memory():
    process = psutil.Process()
    memory_info = process.memory_info()
    return int(memory_info.rss/1024)

# to check for time + execute the sequence alignment
def execute(str_1, str_2):    
    obj = EfficientSequenceAlignment()
    start = time.time()
    n=len(str_2)
    for i in range(0,2):
        opt1.append([0] * (n + 1))
    results = obj.DC_sol(str_1, str_2)    
    end = time.time()
    total_time = (end - start) * 1000

    total_memory = process_memory()

    return [results, total_memory, total_time]

def efficientImplementor(inputFileLocation, outputFileLocation):
    obj = EfficientSequenceAlignment()
    li = obj.input_generate(inputFileLocation)
    
    strings = []
    for i in range(0,len(li)):
        if(li[i].isdigit()==False):
            strings.append(li[i])

    str_1 = strings[0]
    str_2 = strings[1]

    output = execute(str_1, str_2)
    results = output[0]

    f = open(outputFileLocation, 'w')
    f.write(str(int(results[2])) + "\n")
    f.write(results[0][:] + "\n")
    f.write(results[1][:] + "\n")
    f.write(str(float("{:.4f}".format(output[2]))) + "\n") #time
    f.write(str(float(output[1]))) #memory
    f.close()

if __name__ == "__main__":
    efficientImplementor(sys.argv[-2],sys.argv[-1])