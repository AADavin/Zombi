#!/usr/bin/python
# -*- coding: utf-8 -*-

#########################################
##  Author:         Wandrille Duchemin
##  Created:        13-Jan-2017
##  Last modified:  10-Apr-2018
##
##  Decribes 3 classes : RecEvent, ReconciledTree and ReconciledTreeList
##  the ReconciledTree class represent a reconciled gene tree and
##  offers input and output functions for the recPhyloXML format.
##  The ReconciledTreeList class is a container for several instances of ReconciledTree
##  and a facultative species tree
##
##  requires : ete3 ( http://etetoolkit.org/ )
##             xml ( in standard library )
##
##  developped for python2.7
##
#########################################


import ete3
import xml.etree.ElementTree as ET


RECPHYLOTAG = "recPhylo"
RECTREETAG = "recGeneTree"
SPTREETAG = "spTree"


## helper function to get XML lines for a simple XML tree. typically used for the species tree here
def myBasicTreeXMLLinesAux(tree):
    """
    Takes:
        - tree (ete3.TreeNode)
    Returns:
        (list): list of xml lines
    """

    indentChar = "  "

    lines = ["<clade>"]

    lines.append( indentChar + "<name>" + tree.name + "</name>" )

    for c in tree.children:
        tmp = myBasicTreeXMLLinesAux(c)
        for l in tmp:
            lines.append( indentChar + l )

    lines.append("</clade>")

    return lines

def myBasicTreeXMLLines(tree):
    """
    Takes:
        - tree (ete3.TreeNode)
    Returns:
        (list): list of xml lines
    """
    lines = ["<phylogeny>"]
    indentChar = "  "
    tmp = myBasicTreeXMLLinesAux(tree)
    for l in tmp:
            lines.append( indentChar + l )

    lines.append("</phylogeny>")

    return lines


EVENTTAGCORRESPONDANCE = {    "D" : "duplication",
                        "S" : "speciation",
                        "C" : "leaf",
                        "L":"loss",
                        "Bo": "bifurcationOut",
                        "bro": "branchingOut",
                        "Tb": "transferBack",
                        "SL": "speciationLoss",
                        "broL": "branchingOutLoss"
                        }

class RecEvent:
    def __init__(self, eventCode , species , ts = None , additionnalInfo = {}):
        """
        Takes:
            - eventCode (str) : a code indicating the recEvent event
            - species (~) : a identifier for the specie the event takes place in
            - ts (int or None) [default= None] : the time slice the events happens at, if applicable
            - additionnalInfo (dict) [default= {}] : keys are expected to be some property tag and values the associated information
        """
        self.eventCode = eventCode
        self.species = species
        self.timeSlice = ts
        self.additionnalInfo = additionnalInfo.copy()

    def __str__(self):
        eventName = self.eventCode
#        if EVENTTAGCORRESPONDANCE.has_key(self.eventCode):
        if self.eventCode in EVENTTAGCORRESPONDANCE:
            eventName = EVENTTAGCORRESPONDANCE[self.eventCode]

        L = [ str(self.eventCode), "spe=" + str(self.species) ]
        if not self.timeSlice is None:
            L.append("ts=" + str(self.timeSlice))
        for k,v in self.additionnalInfo.items():
            L.append(str(k) + "=" + str(v) )
        return " ".join(L)

    def nwkstr(self):
        """ tmp simplistic version """
        s = str(self.species)
        s += "."

        s += str(self.eventCode)
        return s


    def makeRecXMLstr(self , speciesNames):

        if self.eventCode == "N":
            return ""

        eventName = self.eventCode
#        if EVENTTAGCORRESPONDANCE.has_key(self.eventCode):
        if self.eventCode in EVENTTAGCORRESPONDANCE:
            eventName = EVENTTAGCORRESPONDANCE[self.eventCode]

        spe = str(self.species)
        if spe in speciesNames.keys():
            spe = speciesNames[spe]

        S = "<" + str(eventName)

        if self.eventCode != "Bo":
            S += " "
            if self.eventCode == "Tb":
                S += "destinationSpecies="
            else:
                S += "speciesLocation="
            S += '"' + str(spe) + '"'

            if not self.timeSlice is None:
                S += " ts=" + '"' + str(self.timeSlice) + '"'

        propertyName = "confidence"
#        if self.additionnalInfo.has_key(propertyName):
        if propertyName in self.additionnalInfo:
            S += " " + propertyName + "=" + '"' + self.additionnalInfo[propertyName] + '"'

        if self.eventCode == "C":
            propertyName = "geneName"
#            if self.additionnalInfo.has_key(propertyName):
            if propertyName in self.additionnalInfo:
                S += " " + propertyName + "=" + '"' + self.additionnalInfo[propertyName] + '"'


        S += ">"
        S += "</" + eventName + ">"

        #print self.eventCode , "->", S
        return S


class ReconciledTree(ete3.TreeNode):
    def __init__(self):
        ete3.TreeNode.__init__(self)

        self.name = ""
        self.eventRecs = []

    def setName(self, name):
        """ NB : name can be any object with a __str__() fc (like an int for instance) """
        self.name = name

    def getEvents(self):
        return self.eventRecs

    def getEvent(self,i):
        return self.eventRecs[i]

    def popEvent(self , i):
        return self.eventRecs.pop(i)


    def addEvent(self , e,append = True):
        """
            Takes:
                - e (RecEvent) : the reconciliation event to add
                - append (bool) [default=True] : if True adds the event at the end of the list else, adds it at the beginning
        """
        if append:
            self.eventRecs.append(e)
        else:
            self.eventRecs.insert(0,e)
        return

    def getTreeStrAux(self):

        L = []
        L.append("name : " + str(self.name) )
        L.append("events :")
        for e in self.eventRecs:
            L.append("  " + str(e))
        ChL = []
        for c in self.get_children():
            ChL += ["  " + s for s in c.getTreeStrAux()]
        if len(ChL)>0:
            L.append("children :")
            L += ChL
        return L



    def getTreeStr(self):
        return "\n".join(self.getTreeStrAux())

    def getTreeNewickAux(self, sep="|", topoOnly = False):
        s = ""

        ChL = []
        for c in self.get_children():
            ChL.append(c.getTreeNewickAux(sep,topoOnly))
        if len(ChL)>0:
            s += "("
            s += ",".join(ChL)
            s += ")"

        s += str(self.name)
        if not topoOnly:
            s += sep
            s += sep.join( [e.nwkstr() for e in self.eventRecs] )

        return s

    def getTreeNewick(self, sep="|" , topoOnly = False):
        return self.getTreeNewickAux(sep, topoOnly) + ";"

    def getTreeRecPhyloXMLAux(self , speciesNames ={}, topoOnly = False):

        L = []
        L.append("<clade>")
        L.append("  <name>" + str(self.name) + "</name>" )
        if not topoOnly:
            L.append("  <eventsRec>")
            for e in self.eventRecs:
                s = e.makeRecXMLstr(speciesNames)
                if s == "":
                    continue
                L.append("    " +  s)
            L.append("  </eventsRec>")
        ChL = []
        for c in self.get_children():
            ChL += ["  " + s for s in c.getTreeRecPhyloXMLAux(speciesNames, topoOnly)]
        if len(ChL)>0:
            L += ChL
        L.append("</clade>")
        return L

    def getTreeRecPhyloXML(self , speciesNames ={}, topoOnly = False):
        Lines = self.getTreeRecPhyloXMLLines( speciesNames, topoOnly)
        return "\n".join(Lines)

    def getTreeRecPhyloXMLLines(self , speciesNames ={}, topoOnly = False):
        Lines = ["<recGeneTree>"]
        Lines.append("  <phylogeny rooted=\"true\">")
        tmp = self.getTreeRecPhyloXMLAux( speciesNames , topoOnly )
        for l in tmp:
            Lines.append( "    " + l )
        Lines.append("  </phylogeny>")
        Lines.append( "</recGeneTree>" )
        return Lines

    def countEvents(self):
        """
        Returns:
            (dict) : keys are recPhyloXML event tags, values are the number of times these events occur in the tree
        """
        devent = {}
        for e in self.eventRecs:
            code  = e.eventCode
#            if not devent.has_key(code):
            if not code in devent:
                devent[code] = 0
            devent[code] += 1

        for c in self.get_children():
            tmp = c.countEvents()
            for k in tmp.keys():
#                if not devent.has_key(k):
                if not k in devent:
                    devent[k] = 0
                devent[k] += tmp[k]

        return devent

    def getEventsSummary(self , speciesTree, includeTransferReception  = True , includeTransferDeparture = False , speciesIdFeature = "name"):
        """
        *recursive function*
        Takes:
             - speciesTree (ete3.Tree) : the species tree used for the reconciliation, necessary to assign a species to loss events.
             - includeTransferReception [default = True]  : Whether or not to includes events of  TransferReception (transferBack tag) in the counts.
             - includeTransferDeparture [default = False] : Whether or not to includes events of  TransferDeparture (branchingOut tag) in the counts.
             - speciesIdFeature (str) [default = "name"] : the feature to use as Id in the species tree (by default the name is used)
        Returns:
            (dict):
                    keys are events type among : "duplication" , "loss" , "transferReception" , "transferDeparture"
                    values are  lists of species id
        """

        EventsSummary = { "duplication" : [],
                          "loss" : [] }

        if includeTransferReception:
            EventsSummary["transferReception"] = []

        if includeTransferDeparture:
            EventsSummary["transferDeparture"] = []


        for i,e in enumerate(self.eventRecs):

            evtCode = e.eventCode

            species = e.species

            report = False

            reportTD = False
            reportTR = False

            if ( evtCode in ( EVENTTAGCORRESPONDANCE["SL"] , EVENTTAGCORRESPONDANCE["broL"] ) ) or (evtCode in ( "SL", "broL" )):
                species = self.getLostSpecies( i , speciesTree, speciesIdFeature)
                evtCode = "loss"
                report = True

            elif evtCode == EVENTTAGCORRESPONDANCE["D"] or evtCode =="D":
                evtCode = "duplication"
                report = True

            elif evtCode == EVENTTAGCORRESPONDANCE["L"] or evtCode =="L":
                evtCode = "loss"
                report = True

            if report: ##reporting duplication or loss
                EventsSummary[evtCode].append(species)


            if includeTransferReception and ( evtCode == EVENTTAGCORRESPONDANCE["Tb"] or evtCode == "Tb" ) :
                evtCode = "transferReception"
                EventsSummary[evtCode].append(species)

            elif includeTransferDeparture and ( (evtCode in [ EVENTTAGCORRESPONDANCE["bro"] , EVENTTAGCORRESPONDANCE["broL"] ] ) or ( evtCode in [ "bro", "broL" ]) ):
                evtCode = "transferDeparture"
                EventsSummary[evtCode].append(species)


        for c in self.get_children():

            tmp = c.getEventsSummary(speciesTree , includeTransferReception , includeTransferDeparture , speciesIdFeature)

            for k,v in tmp.items():
                EventsSummary[k] += v

        return EventsSummary




    def sameSpeciesAsParent(self , parent = None):
        """ returns True if the first event of the node has the same species as the last event of its parent , False otherwise (and if self is the root)
            if the parent is given, is it used, otherwise we look for it in the structure
        """

        if parent is None:
            if self.is_root():
                return False
            parent = self.up

        lastParentSp = parent.getEvents()[-1].species

        firstSp = self.getEvents()[0].species

        return firstSp == lastParentSp


    def getLostSpecies( self, evtIndex , speciesTree, speciesIdFeature = "name"):
        """
        given the index of an event of *Loss (speciationLoss for instance) in this nodes,
        this function returns the id of the species where the loss occured
        (this function is useful because speciationLoss event references the species of the speciation rather than the species of the loss)
        Takes:
            - evtIndex (int) : index of the loss event whose lost species we want to know
            - speciesTree (ete3.Tree) : the species tree used for the reconciliation
            - speciesIdFeature (str) [default = "name"] : the feature to use as Id in the species tree (by default the name is used)
        Returns:
            (str) : id of the species where the loss occured
            or
            None : if the indicated event does not correspond to a loss
        """


        evtCode = self.getEvent(evtIndex).eventCode

        if not evtCode in [ "SL","broL",EVENTTAGCORRESPONDANCE["SL"] , EVENTTAGCORRESPONDANCE["broL"] ]:
            return None

        if evtCode in ('broL', EVENTTAGCORRESPONDANCE["broL"]):
            return self.getEvent(evtIndex).species ## in branchingOutLoss, the species of the lost lineage is the same as the one of the branchingOut

        ## We know this is a speciationLoss event.
        ## Note that speciationLoss events are NEVER the last event in eventsRec (as they are neither a bifurcation nor a leaf event).
        ## so we should be able to safely ask for the event after the current one.

        lostSpeciesSister = self.getEvent( evtIndex + 1 ).species

        Speciesnode = speciesTree.search_nodes(**{speciesIdFeature: lostSpeciesSister })

        if len(Speciesnode) != 1:
            raise Error("error:",len(Speciesnode),"with Id",lostSpeciesSister ,"(1 expected).")

        lostSpeciesNode = Speciesnode[0].get_sisters()[0]

        lostSpeciesId = getattr(lostSpeciesNode, speciesIdFeature)


        return lostSpeciesId



class ReconciledTreeList:
    """
    This object represents a group of reconciled tree.
    Usually they would be reconciled with the same species tree, which can also be added to this object.
    Atributes:
        - self.spTree   : the species tree these trees are reconciled with (or None if no species tree is specified)
        - self.recTrees : the reconciled trees
    """
    def __init__(self, spTree = None , recTrees = []):
        """
        Takes:
            spTree (ete3.Tree) [ default = None ] : facultative species tree
            recTrees (list) [ default = [] ] : list of ReconciledTree instance (see above)
        """

        self.spTree = spTree
        self.recTrees = recTrees[:]


    def setSpTree(self, ST):
        """
        Simply sets a trees as the object species tree.
        Takes:
            - ST (ete3.Tree) : a species tree
        """
        self.spTree = ST



    def append(self, RT):
        """
        Appends a reconciled tree to the object.
        Takes:
            - RT (ReconciledTree) : a reconciled tree
        """
        self.recTrees.append(RT)

    def __getitem__(self, i ):
        """
        returns a reconciled tree from the object.
        Takes:
            - i (int) : index of the desired reconciled tree
        Returns:
            (ReconciledTree) : the reconciled tree at the desired index
            OR IndexError if the index is invalid (ie. too high)
        """

        if i >= len(self.recTrees):
            raise IndexError('Index out of range. There are no reconciled tree with index ' + str(i) + '.')

        return self.recTrees[i]


    def __len__(self):
        """
        Returns:
            (int) : number of ReconciledTree in this instance
        """

        return len(self.recTrees)

    def hasSpTree(self):
        """
        Returns:
            (bool) : True if there this instance has a species tree (ie. self.spTree is not None), False otherwise
        """
        return not self.spTree is None

    def getRecPhyloXMLLines(self):
        """
        Returns:
            (list) : list of lines of the recPhyloXML representation of this object
                     (NB : the lines do not have a '\n' at their end.)
        """
        lines = []
        lines.append("<" + RECPHYLOTAG + ">"  )

        offset = 1
        offsetChar = "  "

        if self.hasSpTree():
            lines.append( offsetChar*offset + "<" + SPTREETAG + ">"  )
            offset += 1
            spLines = myBasicTreeXMLLines(self.spTree)
            for l in spLines:
                lines.append( offsetChar*offset + l )
            offset -= 1
            lines.append( offsetChar*offset + "</" + SPTREETAG + ">"  )


        for RT in self.recTrees:
            recLines = RT.getTreeRecPhyloXMLLines()
            for l in recLines:
                lines.append( offsetChar*offset + l )


        lines.append("</" + RECPHYLOTAG + ">" )

        return lines

    def getEventsSummary(self , includeTransferReception  = True , includeTransferDeparture = False , speciesIdFeature = "name" , indexBySpecies=False):
        """
        Retrieve an event summary over all the trees in the object
        !!only works if there is a species tree assigned to the object!!
        Takes:
             - includeTransferReception [default = True]  : whether or not to includes events of  TransferReception (transferBack tag) in the counts.
             - includeTransferDeparture [default = False] : whether or not to includes events of  TransferDeparture (branchingOut tag) in the counts.
             - speciesIdFeature (str) [default = "name"] : the feature to use as Id in the species tree (by default the name is used)
             - indexBySpecies (str) [default = False] : if True, the returned dictionnary will have species as keys and event counts as values.
        Returns:
            (dict):
                    keys are events type among : "duplication" , "loss" , "transferReception" , "transferDeparture"
                    values are  lists of species id
                   OR, if indexBySpecies=True:
                       keys are species id
                       values are dict with keys among "duplication" , "loss" , "transferReception" , "transferDeparture"
                                            and values as counts of the events in each species
        """

        if not self.hasSpTree():
            raise Error("error : can't get an events summary when no species tree has been assigned.")

        EventsSummary = { "duplication" : [],
                          "loss" : [] }

        if includeTransferReception:
            EventsSummary["transferReception"] = []

        if includeTransferDeparture:
            EventsSummary["transferDeparture"] = []

        for RT in self.recTrees:
            tmp = RT.getEventsSummary(self.spTree , includeTransferReception , includeTransferDeparture , speciesIdFeature)



            for k,v in tmp.items():
                EventsSummary[k] += v




        if indexBySpecies:

            tmp = {}

            for n in self.spTree.traverse():

                tmp[ getattr(n,speciesIdFeature) ] = {}



            for e in EventsSummary.keys():

                for sp in tmp.keys():
                    tmp[sp][e] = 0

                for sp in EventsSummary[e]:
                    tmp[sp][e] += 1

            EventsSummary = tmp


        return EventsSummary