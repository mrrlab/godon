package cmodel

import (
	"errors"
	"fmt"
	"os"

	"bitbucket.org/Davydov/godon/bio"
	"bitbucket.org/Davydov/godon/codon"
	"bitbucket.org/Davydov/godon/tree"
)

// Data stores data for the codon model
type Data struct {
	// cSeqs is the codon alignment
	cSeqs codon.Sequences
	// tree is the tree
	Tree *tree.Tree
	// cFreq is codon frequencies
	cFreq codon.Frequency
	// if the tree was rooted
	root bool
	// old root node id
	rootID int
}

// NewData creates a new Data.
func NewData(gCodeID int, alignmentFileName string, treeFileName string,
	cFreq string) (*Data, error) {
	data := &Data{}

	gcode, ok := bio.GeneticCodes[gCodeID]
	if !ok {
		return nil, fmt.Errorf("couldn't load genetic code with id=%d", gCodeID)
	}
	log.Infof("Genetic code: %d, \"%s\"", gcode.ID, gcode.Name)

	var fastaFile *os.File
	fastaFile, err := os.Open(alignmentFileName)
	if err != nil {
		return nil, err
	}
	defer fastaFile.Close()

	ali, err := bio.ParseFasta(fastaFile)
	if err != nil {
		return nil, err
	}

	data.cSeqs, err = codon.ToCodonSequences(ali, gcode)
	if err != nil {
		return nil, err
	}

	if data.cSeqs.Length() == 0 {
		return nil, errors.New("Zero length alignment")
	}
	log.Infof("Read alignment of %d codons, %d fixed positions, %d ambiguous positions", data.cSeqs.Length(), data.cSeqs.NFixed(), data.cSeqs.NAmbiguous())

	treeFile, err := os.Open(treeFileName)
	if err != nil {
		return nil, err
	}
	defer treeFile.Close()

	data.Tree, err = tree.ParseNewick(treeFile)
	if err != nil {
		return nil, err
	}

	switch cFreq {
	case "F0":
		log.Info("F0 frequency")
		data.cFreq = codon.F0(data.cSeqs)
	case "F3X4":
		log.Info("F3X4 frequency")
		data.cFreq = codon.F3X4(data.cSeqs)
	default:
		return nil, errors.New("Unknow codon freuquency specification")
	}

	return data, nil
}

// SetForegroundBranch sets the foreground branch in a tree.
func (data *Data) SetForegroundBranch(fgBranch int) {
	for _, node := range data.Tree.NodeIDArray() {
		if node == nil {
			continue
		}
		if node.ID == fgBranch {
			node.Class = 1
		} else {
			node.Class = 0
		}
	}
}

// GetNClass1 returns number of nodes with class=1.
func (data *Data) GetNClass1() (res int) {
	for range data.Tree.ClassNodes(1) {
		res++
	}
	return
}

// SetCodonFreqFromFile sets codon frequency from file.
func (data *Data) SetCodonFreqFromFile(filename string) error {
	cFreqFile, err := os.Open(filename)
	if err != nil {
		return err
	}
	data.cFreq, err = codon.ReadFrequency(cFreqFile, data.cFreq.GCode)
	return err
}

// Copy creates a copy (only new tree is created).
func (data *Data) Copy() *Data {
	return &Data{
		cSeqs:  data.cSeqs,
		Tree:   data.Tree.Copy(),
		cFreq:  data.cFreq,
		root:   data.root,
		rootID: data.rootID,
	}
}

// Unroots the tree if needed.
func (data *Data) Unroot() error {
	var err error
	if data.Tree.IsRooted() {
		log.Warning("Tree is rooted. Will unroot.")
		data.root = true
		data.rootID, err = data.Tree.Unroot()
		if err != nil {
			return fmt.Errorf("Error unrooting tree: %v", err)
		}
	}

	log.Infof("intree_unroot=%s", data.Tree)
	log.Debugf("brtree_unroot=%s", data.Tree.BrString())
	log.Debug(data.Tree.FullString())
	return nil
}

// Root roots the tree back in case it was unrooted.
func (data *Data) Root() error {
	if data.root {
		log.Infof("unrooted_outtree=%s", data.Tree)
		err := data.Tree.Root(data.rootID)
		if err != nil {
			return fmt.Errorf("Error rooting tree: %v", err)
		}
	}
	return nil
}
