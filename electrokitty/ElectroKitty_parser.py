# -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 10:41:11 2024

@author: ozbejv
"""

class electrokitty_parser:
    """
    Class containing the mechanism parser
    
    The parser uses the mechanism string to generate lists used by the simulator
    
    The final output is a list countaining everything the simulator needs
    
    A script of functions for parsing a mechanism string
    Zero library dependancies requierd, for Python 3.x 
    Mechanistic steps must be separated by \n, the species in forward and backward reaction are
    separated by +, with = stating that the reaction is in equilibrium, and - saying 
    the reaction is one-directional
    Besides the characters, \n, +, =, -, any character or sequence of characters can be used,
    With the caveat that all spaces are discraded
    """
    def __init__(self, string):
        self.string=string
        self.mechanism_list=None
    
    ################# Functions for parsing mechanism
    def update_species_lists(self, line, bulk_spec, ads_spec, string):
        """
        For a step in the mechanism sorts the species into either adsorbed or dissolved
        and updates the appropriate lists
        """
        for reac_fb in line.split(string):
            for s in reac_fb.split("+"):
                if "*" == s[-1]:
                    if s not in ads_spec:
                        ads_spec.append(s)
                else:
                    if s not in bulk_spec:
                        bulk_spec.append(s)
        return bulk_spec, ads_spec
    
    def find_string_index(self, string, lis):
        """
        given a string and list finds the index where the sring is located in the list
        """
        i=0
        for element in lis:
            if element == string:
                return i
            else:
                i+=1
    
    def get_index(self, string, master_list):
        """
        wrapper function for previous function that gives the correct index given that he species is either adsorbed or not
        """
        if string in master_list[0]:
            return self.find_string_index(string, master_list[0])
        elif string in master_list[1]:
            return self.find_string_index(string, master_list[1])+len(master_list[0])
        else:
            print("Mechanism error: Cannot find species in species list")
            return False
    
    def index_list(self, line, species_list):
        """
        for each section of the mechanistic step finds the index of the species
        """
        f_index=[]
        
        for spec in line.split("+"):
            f_index.append(self.get_index(spec, species_list))
            
        return f_index
    
    def update_index_lists(self, ads_index, bulk_index , ec_index, f,b, line,
                           reaction_types_ads, reaction_types_ec, reaction_types_bulk, Type, len_ads, reaction_index, indicator, num_el):
        
        """
        given the mechanism step updates the lists of the species, depending on the reaction of the step
        Updates the list of species as well as the indexes for the reaction
        """
        
        if "*" in line and line[:2]=="C:":
            ads_index.append([f,b])
            reaction_types_ads.append(Type)
            reaction_index[0].append(indicator)
        elif line[0]=="E":
            elelectron_number=int(line[2])
            num_el.append(elelectron_number)
            ec_index.append([f,b])
            reaction_types_ec.append(Type)
            reaction_index[2].append(indicator)
        else:
            for i in range(len(f)):
                f[i]-=len_ads
            for i in range(len(b)):
                b[i]-=len_ads
            bulk_index.append([f,b])
            reaction_types_bulk.append(Type)
            reaction_index[1].append(indicator)
        return ads_index, bulk_index, ec_index, reaction_types_ads, reaction_types_ec, reaction_types_bulk, reaction_index, num_el
    
    def Parse_mechanism(self):
        """
        # The main function for analysing the mechanism string.
        # # The function takes the string, with each line (separated by \n),
        # each line must be declared as either a C or E mechanism, followed by :,
        # * at the end of the species, assignes the species as adsorbed,
        # and gives out a list of species, first list beiing the adsorbed, second the soulution,
        # and three lists of indecies to connect the species via differential equetions.
        # First are the adsobed species list, second is the bulk and finally the ec connections
        
        # < backaward step 1    > forward step  2
        """
        
        
        string=self.string
        string=string.replace(" ","")
        a=string.split("\n")
    
        bulk_spec=[]
        ads_spec=[]
    
        bulk_index=[]
        ads_index=[]
        ec_index=[]
        
        reaction_types_ads=[]
        reaction_types_ec=[]
        reaction_types_bulk=[]
        
        num_el=[]
        
        for line in a:
            
            ind=self.find_string_index(":", line)
            if "=" in line:
                bulk_spec, ads_spec = self.update_species_lists(line[ind+1:], bulk_spec, ads_spec, "=")
            elif "<" in line:
                bulk_spec, ads_spec = self.update_species_lists(line[ind+1:], bulk_spec, ads_spec, "<")
            elif ">" in line:
                bulk_spec, ads_spec = self.update_species_lists(line[ind+1:], bulk_spec, ads_spec, ">")
            else:
                print("Mechanism Error: Wrong mechanism separator")
    
        species_list=[ads_spec,bulk_spec]
        
        reaction_index=[[],[],[]]
        for i in range(len(a)):
            line=a[i]
            ind=self.find_string_index(":", line)
            if "=" in line:
                f,b=line[ind+1:].split("=")
                f_ind=self.index_list(f, species_list)
                b_ind=self.index_list(b, species_list)
                ads_index, bulk_index, ec_index, reaction_types_ads,reaction_types_ec,reaction_types_bulk,reaction_index,num_el = self.update_index_lists(
                    ads_index, bulk_index, ec_index, f_ind, b_ind, line,
                    reaction_types_ads,
                    reaction_types_ec,reaction_types_bulk,0, len(species_list[0]), reaction_index, i, num_el)
                
            elif "<" in line:
                f,b=line[ind+1:].split("<")
                f_ind=self.index_list(f, species_list)
                b_ind=self.index_list(b, species_list)
                ads_index, bulk_index, ec_index, reaction_types_ads,reaction_types_ec,reaction_types_bulk,reaction_index,num_el = self.update_index_lists(
                    ads_index, bulk_index, ec_index, f_ind, b_ind, line,
                    reaction_types_ads,
                    reaction_types_ec,reaction_types_bulk,1, len(species_list[0]), reaction_index, i, num_el)
            elif ">" in line:
                f,b=line[ind+1:].split(">")
                f_ind=self.index_list(f, species_list)
                b_ind=self.index_list(b, species_list)
                ads_index, bulk_index, ec_index, reaction_types_ads,reaction_types_ec,reaction_types_bulk,reaction_index,num_el = self.update_index_lists(
                    ads_index, bulk_index, ec_index, f_ind, b_ind, line,
                    reaction_types_ads,
                    reaction_types_ec,reaction_types_bulk,2, len(species_list[0]), reaction_index, i, num_el)
            else:
                print("Mechanism error: cannot index mechanism")
        self.mechanism_list=[species_list, [ads_index, bulk_index, ec_index], [reaction_types_ads, reaction_types_bulk, reaction_types_ec], reaction_index, num_el]
        return self.mechanism_list