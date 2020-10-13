from __future__ import print_function
import subprocess
import re
import sys
from collections import OrderedDict

baseFolder = sys.argv[1]

print('/**')
print('\page parameters Ptree-based Parameters')
print('Click on values below to expand parameter tree and see details.  The "key" value of the details can be used with the ptree::put() or ptree::get() functions to access the value in a boost::property_tree::ptree instance.')
print(' ')
print('\htmlonly')

f = open(baseFolder+'/documentation/doxFiles/listHeader.txt','r')
headerLines = f.read()
print(headerLines)
f.close()

print('\endhtmlonly')

p = subprocess.Popen(["find",baseFolder+"/modules","-type","f", "-name","*.cpp", "-exec", "grep","-ni","\.get(\|ReadAndLogParameter","{}","+"],stdout=subprocess.PIPE)

out,err = p.communicate()
if(out==b''):
    filenames = []
    allData = []

else:

    out = str(out).split('\n')
    filenames = set([str(x.split('.cpp:',1)[0]+".cpp") for x in out])
    filenames.discard('.cpp')

    allData = []

    fileList = []
    lineNumList = []
    varList = []
    defaultList = []

    # now that we have all the files and use ptree.get.  Open them up and read the ptree values
    for file in filenames:

        if(('test' not in file)&(len(file)>4)):
            try:
                f = open(file,'r')
                allLines = f.read()
                f.close()

                for m in re.finditer('\.get\("|ReadAndLogParameter',allLines):
                    startInd = m.start()
                    endInd=startInd
                    openBrackets = 0
                    closedBrackets = 0
                    while((openBrackets!=closedBrackets)|(openBrackets==0)):
                        endInd += 1
                        if(allLines[endInd]=='('):
                            openBrackets += 1
                        elif(allLines[endInd]==')'):
                            closedBrackets += 1

                    lineNum = allLines[0:startInd].count('\n')+1

                    tempInd = file.find('MUQ/modules')

                    fileList.append(file[tempInd:])
                    lineNumList.append(lineNum)

                    unpacked = allLines[startInd:endInd].replace("\n","").split(',')
                    if(len(unpacked)==2):
                        [variable, default] = allLines[startInd+5:endInd].replace("\n","").split(',')
                        variable = variable.strip().strip('"')
                        default = default.strip()
                        default.replace(" ","")

                    else:
                        variable = unpacked[1].strip().strip('"')
                        default = unpacked[2].strip().replace(" ","")

                    varList.append(variable.strip().strip('"'))
                    defaultList.append(default)
            except:
                print('Could not open %s' + file)

    allData = list(zip(varList,defaultList,lineNumList, fileList))
    allData.sort()

graph = OrderedDict()
data = dict()
isRoot = dict()

for x in allData:
    allNodes = x[0].split('.')
    data[x[0]] = x[1:]

    for j in range(0,len(allNodes)):
        key = '.'.join(allNodes[0:j])
        child = '.'.join(allNodes[0:j+1])

        if(key not in graph.keys()):
            graph[key] = []
        if(child not in graph.keys()):
            graph[child] = []

        graph[key].append(child)
    graph[x[0]] = []

def PrintChild(node, graph, data, isLast):
  if(isLast):
    print('<li class="lastChild">', node.split('.')[-1])
    return
  else:
    print("<li>", node.split('.')[-1])
  if(len(graph[node])==0):
    print("<ul>")
    print("<li> Default Value = ", data[node][0], "</li>")
    print('<li> key = "' + node + '"</li>')
    nameSplit = data[node][2].split('/')[-1].split('.')
    fileLink = nameSplit[0]+'_8'+nameSplit[1]+'_source.html'
    lineLink = fileLink + '#l' + str(data[node][1]).zfill(5)
    print('<li class="lastChild"> Extracted from line <a href=' + lineLink +'>'+ str(data[node][1]) + "</a> of <a href=" + fileLink+'>'+data[node][2]+'</a>',"</li>")
    print("</ul>")
  else:
    print('<ul class="collapsibleList">')

    temp = sorted(list(set(graph[node])))
    for child in range(len(temp)-1):
      PrintChild(temp[child],graph,data,0)
    PrintChild(temp[len(temp)-1],graph,data,1)

    print("</ul>")
  print("</li>")


print('<ul class="collapsibleList">')
for node in graph:
    if((node.count('.')==0)&(len(node)>0)):
        PrintChild(node,graph,data,0)
print('</ul>')
print(' ')
print('*/')
print(' ')
