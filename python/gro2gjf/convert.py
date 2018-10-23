#!/usr/bin/env python
# -*- coding:utf-8 -*-
#Author:firestige

import re, os
from pathlib import Path
from decimal import Decimal

class Element(object):

    def __init__(self, *agrs):
        self.e = agrs[0]
        self.index = int(agrs[1])
        self.x = Decimal(agrs[2])*10
        self.y = Decimal(agrs[3])*10
        self.z = Decimal(agrs[4])*10

    def export2Str(self, limit):
        if limit:
            return "{:>3}{:>6}{:>25.14f}{:>25.14f}{:>25.14f}\n".format(self.e, "-1", self.x, self.y, self.z)
        else:
            return "{:>3}{:>6}{:>25.14f}{:>25.14f}{:>25.14f}\n".format(self.e, "0", self.x, self.y, self.z)

    def getSqrtDistance(self, centre):
        return (self.x-centre.x)*(self.x-centre.x) + (self.y-centre.y)*(self.y-centre.y) + (self.z-centre.z)*(self.z-centre.z)

class GroHandler(object):

    def __init__(self, inputPath, outputPath, centre_index, distance, start_index, protein_atom_nums, limit_atom_index):
        #输入文件
        inputPath = Path(inputPath)
        if inputPath.is_file():
            self.inputFile = open(inputPath,"r", encoding="utf-8")
        else:
            raise FileNotFoundError
        #输出文件
        outputPath = Path(outputPath)
        if outputPath.is_file():
            os.remove(outputPath)
        self.outputFile = outputPath
        #中心碳原子序号
        self.centre_index = centre_index
        #标准距离
        self.distance = distance
        #坐标起始行位置
        self.start_index = start_index
        #蛋白质原子数
        self.protein_atom_nums = protein_atom_nums
        #被限制位置的原子序号
        self.limit_atom_index = limit_atom_index
 
    def write2file(self, content):
            p = Path(self.outputFile)
            if p.is_file():
                with p.open(mode='a',encoding='utf-8') as f:
                    f.write(content)
            else:
                with p.open(mode='w',encoding='utf-8') as f:
                    f.write(content)

    def execute(self):        
        #初始化中心元素
        centre = None
        #计数器归零
        index = 0
        #前进至起始行
        while index < self.start_index:
            next(self.inputFile)
            index += 1
        #计数器归零    
        index = 0
        #录入氨基酸小分子
        while index < self.protein_atom_nums:
            element = Element(*re.split(r"\s+", next(self.inputFile).strip())[1:6])
            if element.index == self.centre_index:
                centre = element
                self.write2file(element.export2Str(True))
            elif element.index in self.limit_atom_index:
                self.write2file(element.export2Str(True))
            else:
                self.write2file(element.export2Str(False))
            index += 1
        #录入水
        try:
            while True:
                line = next(self.inputFile)
                if len(re.split(r"\s+", line.strip())) < 9:
                    raise StopIteration
                element = Element(*re.split(r"\s+", line.strip())[1:6])
                if centre != None and len(element.e) > 1 and "O" in element.e and element.getSqrtDistance(centre) < self.distance:
                    eO = element.export2Str(True)
                    eH1 = Element(*re.split(r"\s+", next(self.inputFile).strip())[1:6]).export2Str(True)
                    eH2 = Element(*re.split(r"\s+", next(self.inputFile).strip())[1:6]).export2Str(True)
                    next(self.inputFile)
                    self.write2file(eO+eH1+eH2)
        except StopIteration:
            print("文件录入完毕")
        finally:
            self.inputFile.close()

if __name__ == "__main__":
    
    #单独调用时按照下述方法调用
    '''
    setting = {
        "inputPath":"npt9.gro",#输入文件，type=str
        "outputPath":"temp.xyz",#输出文件，type=str
        "centre_index":5,#被固定的原子序号，type=int
        "distance":25,#选择半径，type=int
        "start_index":2,#忽略行数，type=int
        "protein_atom_nums":12#蛋白质小分子原子数，type=int
    }
    handler = GroHandler(**setting)
    handler.execute()
    '''
    
    #遍历输入文件夹下的gro文件，并在输出文件夹下创建同名gjf文件
    inputDirPath = "./input/"
    outputDirPath = "./output"
    for inputfile in Path(inputDirPath).rglob('*.gro'):
        setting = {
            "inputPath":inputfile,#输入文件路径，type=str
            "outputPath":outputDirPath+"/"+inputfile.stem+".gjf",#输出文件路径，type=str
            "centre_index":5,#被固定的原子序号，type=int
            "distance":25,#选择半径，type=int
            "start_index":2,#忽略行数，type=int
            "protein_atom_nums":12,#蛋白质小分子原子数，type=int
            "limit_atom_index":[5]#被限定位置的原子序号，type=list
        }
        handler = GroHandler(**setting)
        handler.execute()
