#!/usr/bin/env python
# coding: utf-8

# In[14]:


from multiprocessing import Pool

import pandas as pd
from bs4 import BeautifulSoup
from selenium import webdriver
from selenium.webdriver.support.ui import Select
from selenium.common.exceptions import NoSuchElementException


# In[15]:


class Parser():
    
    # Open a Chrome broswer
    def __init__(self, startIndex, process):
        self.browser = webdriver.Chrome()
        self.process = process
        self.numHit = 100
        self.numTotalData = 14479
        self.numPage = self.numTotalData // self.numHit + 1
        self.numData = startIndex
        self.index = self.numData % self.numHit
        self.pIndex = self.numData // self.numHit
        self.dict_nmr_data = dict()
        
        try:
            self.data = pd.read_pickle("nmrDB"+str(self.process)+".pkl")
        except FileNotFoundError:
            self.data = pd.DataFrame()
        
        try:
            self.exceptions = pd.read_pickle("exceptions"+str(self.process)+".pkl")
        except FileNotFoundError:
            self.exceptions = pd.DataFrame()
        
    # Parse every NMR data in SDBS
    def parseNMRData(self):
        self.navigateToMain()
        self.navigateToList(self.pIndex)
        while self.numData < self.numTotalData:
            self.pIndex = self.numData//self.numHit
            while self.pIndex < self.numPage:                
                if not self.parseList(self.pIndex):
                    self.navigateToList(self.pIndex)
                    break
                if not self.pIndex == self.numPage-1:
                    if self.navigateToNextList():
                        self.navigateToList(self.pIndex)
                        break
                self.pIndex += 1
            
        
    # Return True if the navigated page is asking for an agreement to disclaimer
    def checkDisclaimer(self):
        soup = BeautifulSoup(self.browser.page_source, 'html.parser')
        return "/sdbs/cgi-bin/cre_disclaimer.cgi?REQURL=/sdbs/cgi-bin/direct_frame_top.cgi&amp;amp;REFURL=" in soup.text
    
    # Navigate to the main page of SDBS
    def navigateToMain(self):
        self.browser.get("https://sdbs.db.aist.go.jp/sdbs/cgi-bin/direct_frame_top.cgi")
    
    # Click the agree button of the disclaimer page
    def agreeDisclaimer(self):
        self.browser.switch_to.default_content()
        self.browser.switch_to.frame("Down")
        agreeButton = self.browser.find_element_by_xpath("//input[@type='submit']")
        agreeButton.click()
        
    # Search for the 1H NMR data in the main page of SDBS
    def search1HNMR(self):
        self.browser.switch_to.default_content()
        self.browser.switch_to.frame("Down")
        hitComboBox = Select(self.browser.find_element_by_xpath("/html/body/form/center/table/tbody/tr[2]/td/select[1]"))
        hitComboBox.select_by_value(str(self.numHit))
        self.browser.find_element_by_xpath("/html/body/form/center/table/tbody/tr[1]/td[3]/table/tbody/tr[1]/td/table/tbody/tr[3]/td[1]/input").click()
        self.browser.find_element_by_xpath("/html/body/form/center/table/tbody/tr[2]/td/input[1]").click()
    
    # Parse the NMR data in a pIndex-th list page
    def parseList(self, pIndex):
        self.index = self.numData%self.numHit
        while self.index < self.numHit and self.numData < self.numTotalData:
            if not self.parseMol(self.index):
                return False
        self.index = 0
        return True
    
    # Navigate ot the pIndex-th list page from the disclaimer page
    def navigateToList(self, pIndex):
        self.navigateToMain()
        if self.checkDisclaimer():
            self.agreeDisclaimer()
        self.search1HNMR()
        
        numNav = 0 if pIndex < 30 else pIndex//15-1
        for i in range(numNav):
            self.browser.switch_to.default_content()
            self.browser.switch_to.frame("Down")
            if i == 0:
                navButton = self.browser.find_element_by_xpath("/html/body/center/table/tbody/tr["+str(self.numHit+4)+"]/td/a[29]")
            else:
                navButton = self.browser.find_element_by_xpath("/html/body/center/table/tbody/tr["+str(self.numHit+4)+"]/td/a[30]")
            navButton.click()
            if self.checkDisclaimer():
                return False
        if pIndex > 0 and pIndex < 30:
            listButton = self.browser.find_element_by_xpath("/html/body/center/table/tbody/tr["+str(self.numHit+4)+"]/td/a["+str(pIndex)+"]")
        elif pIndex >= self.numPage - (self.numPage % 15):
            listButton = self.browser.find_element_by_xpath("/html/body/center/table/tbody/tr["+str(self.numHit+4)+"]/td/a["+str(pIndex%15+21)+"]")
        elif pIndex >= 30 and pIndex < self.numPage - (self.numPage % 15):
            if pIndex%15 != 14:
                listButton = self.browser.find_element_by_xpath("/html/body/center/table/tbody/tr["+str(self.numHit+4)+"]/td/a["+str(16+pIndex%15)+"]")
            else:
                listButton = self.browser.find_element_by_xpath("/html/body/center/table/tbody/tr["+str(self.numHit+4)+"]/td/b/font")
        else:
            listButton = self.browser.find_element_by_xpath("/html/body/center/table/tbody/tr["+str(self.numHit+4)+"]/td/b/font")
        listButton.click()
        return self.checkDisclaimer()
    
    # Navigate to the next list page
    def navigateToNextList(self):
        self.browser.switch_to.default_content()
        self.browser.switch_to.frame("Down")
        if self.pIndex >= self.numPage - 15:
            nextListButton = self.browser.find_element_by_xpath("/html/body/center/table/tbody/tr[104]/td/a["+str(self.numPage-self.pIndex+2)+"]")
        elif self.pIndex < 15:
            nextListButton = self.browser.find_element_by_xpath("/html/body/center/table/tbody/tr["+str(self.numHit+4)+"]/td/a["+str(self.pIndex+1)+"]")
        else:
            nextListButton = self.browser.find_element_by_xpath("/html/body/center/table/tbody/tr["+str(self.numHit+4)+"]/td/a[16]")
        nextListButton.click()
        return self.checkDisclaimer()

    def parseMol(self, index):
        
        self.navigateToDetail(index)
        self.parseMainPeak()
        self.navigateToPeak()
        self.parsePeak()
        self.index += 1
        self.numData += 1
        self.navigateToResult()
        return True
                
    # Navigate to the detail page of index-th molecule in the list page
    def navigateToDetail(self, index):
        self.browser.switch_to.default_content()
        self.browser.switch_to.frame("Down")
        detailButton = self.browser.find_element_by_xpath("/html/body/center/table/tbody/tr["+str(3+index)+"]/td[9]/a")
        detailButton.click()
        return self.checkDisclaimer()
    
    def parseMainPeak(self):
        self.dict_nmr_data = dict()
        
        self.browser.switch_to.default_content()
        self.browser.switch_to.frame("Down")
        frame = self.browser.find_element_by_xpath("/html/frameset/frame[2]")
        self.browser.switch_to.frame(frame)
        
        mainPeaks = list()
        try:
            lines = self.browser.find_element_by_xpath("/html/body/table[2]/tbody/tr[2]/td/pre").text.splitlines()
        except:
            lines = self.browser.find_element_by_xpath("/html/body/table[2]/tbody/tr/td[1]/pre").text.splitlines()
        for line in lines:
            elems = line.split()
            if len(elems) == 2:
                try:
                    float(elems[1])
                    mainPeaks.append(elems)
                except ValueError:
                    pass
        
        self.dict_nmr_data['mainPeaks'] = mainPeaks
    
    # Navigate to the peak data page in the detail page
    def navigateToPeak(self):
        self.browser.switch_to.default_content()
        self.browser.switch_to.frame("Down")
        frame = self.browser.find_element_by_xpath("/html/frameset/frame[2]")
        self.browser.switch_to.frame(frame)
        peakButton = self.browser.find_element_by_xpath("/html/body/form/input[5]")
        peakButton.click()
        return self.checkDisclaimer()
    
    # Parse data as a dictionary and append it to the self.data
    def parsePeak(self):        
        soup = BeautifulSoup(self.browser.page_source, 'html.parser')
        
        # Parse peak data
        peakData = soup.select('body > pre')[1].get_text().split()
        peaks = list()
        for i in range(0, len(peakData), 3):
            peaks.append(peakData[i:i+3])
        self.dict_nmr_data['peaks'] = peaks
        # Parse solvent data
        self.dict_nmr_data['solvent'] = soup.select('body > table > tbody > tr > td')[3].text[:-1]
        # Parse InChI code
        for candidate in soup.select('body > table > tbody > tr > td'):
            if candidate.text.strip()[:6] == "InChI=":
                self.dict_nmr_data['inchi'] = candidate.text.strip()
                break
        # Parse molecule name
        try:
            self.dict_nmr_data['name'] = soup.select('body > table > tbody > tr > td > font > b')[1].text[:-1]
        except IndexError:
            self.dict_nmr_data['name'] = "None"
            
        self.data = self.data.append(self.dict_nmr_data, ignore_index=True)
        self.data.to_pickle("nmrDB"+str(self.process)+".pkl")
        print("Parsed "+ str(self.numData) + "th molecule: " + self.dict_nmr_data['name'])
        
    def navigateToResult(self):
        self.browser.switch_to.default_content()
        self.browser.switch_to.frame("Down")
        frame = self.browser.find_element_by_xpath("/html/frameset/frame[1]")
        self.browser.switch_to.frame(frame)
        returnButton = self.browser.find_element_by_xpath("/html/body/a[2]")
        returnButton.click()
        return self.checkDisclaimer()


# In[16]:


def parseNMRDB(info):
    startIndex, endIndex, process = info
    parser = Parser(startIndex, process)
    parser.numData = startIndex
    while parser.numData < endIndex:
        try:
            parser.parseNMRData()
        except NoSuchElementException:
            parser.exceptions.append({"num": parser.numData, "page": parser.pIndex, "index": parser.index}, ignore_index=True)
            parser.exceptions.to_pickle("exceptions"+str(process)+".pkl")
            print("Cannot parse " + str(parser.numData) + "th molecule.")

            parser.index += 1
            parser.numData += 1


# In[17]:

if __name__ == '__main__':
	infos = [[0, 4000, 0], [4000, 8000, 1], [8000, 12000, 2], [12000, 14479, 3]]
	
	pool = Pool(processes=4)
	pool.map(parseNMRDB, infos)


# In[ ]:




