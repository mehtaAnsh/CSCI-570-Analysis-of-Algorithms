import time
import sys
import psutil

class BasicSequenceAlignment:
    def alignSequence(self, str_1, str_2):
        
        l_s1 = len(str_1)
        l_s2 = len(str_2)
        dp=[]
        x = ""
        y = ""
        
        for i in range(l_s2+1):
            li=[]
            for j in range(l_s1+1):
                li.append(0)
            dp.append(li)
            
        alpha = {'AA': 0, 'CC': 0, 'GG': 0, 'TT': 0, 
                'AC': 110, 'CA': 110, 'AG': 48, 'GA': 48, 
                'AT': 94, 'TA': 94, 'CG': 118, 'GC': 118,
                'CT': 48, 'TC': 48, 'GT': 110, 'TG': 110}
        delta = 30

        # initialize values in the DP matrix   
        for i in range(1, l_s1 + 1):
            dp[i][0] = delta*i
        for j in range(0, l_s2 + 1):
            dp[0][j] = delta*j

        for i in range(1, l_s1 + 1):
            for j in range(1, l_s2 + 1):
                string = str_1[i - 1] + str_2[j - 1]
                dp[i][j] = min(dp[i - 1][j - 1] + alpha[string], dp[i][j-1] + delta, dp[i-1][j] + delta)

        i=l_s1 
        j=l_s2

        while i and j:
            st = str_1[i - 1] + str_2[j - 1]
            
            if dp[i][j] == dp[i - 1][j] + delta:
                x = str_1[i - 1] + x
                y = "_" + y
                i -= 1
            elif dp[i][j] == dp[i][j - 1] + delta:
                y = str_2[j - 1] + y
                x = "_" + x
                j -= 1
            elif dp[i][j] == dp[i-1][j-1] + alpha[st]:
                x = str_1[i - 1] + x
                y = str_2[j - 1] + y
                i -= 1
                j -= 1
        
        while i:
            x = str_1[i - 1] + x
            y = "_" + y
            i -= 1

        while j:
            x = "_" + x
            y = str_2[j - 1] + y
            j -= 1

        return x, y, dp[l_s1][l_s2]

    def input_generate(self, inputFileLocation):
        li = []

        with open(inputFileLocation) as f:
            lines = [x.rstrip("\n") for x in f.readlines()]
            for j in range(0, len(lines)):
                if(lines[j].isdigit()==False):
                    strind=j 
                else:
                    k = int(lines[j])
                    strr = lines[strind][0:k+1] + lines[strind] + lines[strind][k+1:len(lines[strind][:])]        
                    lines[strind] = strr
            li.append(lines)

        li = li[0]
        f.close()
        return li

# to check for memory
def process_memory():    
    process = psutil.Process()
    memory_info = process.memory_info()
    return int(memory_info.rss/1024)

# to check for time + execute the sequence alignment
def execute(str_1, str_2):    
    obj = BasicSequenceAlignment()

    start = time.time()
    results = obj.alignSequence(str_1, str_2)    
    end = time.time()
    total_time = (end - start) * 1000

    total_memory = process_memory()

    return [results, total_memory, total_time]

def basicImplementor(inputFileLocation, outputFileLocation):
    
    obj = BasicSequenceAlignment()
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
    f.write(str(float("{:.4f}".format(output[2]))) + "\n") # time
    f.write(str(float(output[1]))) # memory
    f.close()
    
if __name__ == '__main__':
    basicImplementor(sys.argv[-2],sys.argv[-1])