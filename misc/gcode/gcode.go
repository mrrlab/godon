// gcode is a tool to generate files with genetic code in go format
// from asn1 file.
//
// More information is available here:
// - https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
// - ftp://ftp.ncbi.nih.gov/entrez/misc/data/gc.prt
package main

import (
	"bufio"
	"errors"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"
	"unicode"
)

type Mode int

const (
	// Default parsing mode.
	NORMAL Mode = iota
	TABLE
	ASSIGN
	LIST
	ELEMENT
	ELEMENTPAR
	ELEMENTPRECOMMA
	PRECOMMA
	END
)

type GeneticCode struct {
	Name      string
	ShortName string
	Id        int
	Ncbieaa   string
	Sncbieaa  string
}

func (gc GeneticCode) String() string {
	return fmt.Sprintf("<GC: Name=\"%s\", ShortName=\"%s\", Id=%d, A=\"%s\", S=\"%s\">",
		gc.Name, gc.ShortName, gc.Id, gc.Ncbieaa, gc.Sncbieaa)
}

func (gc GeneticCode) GoString() string {
	return fmt.Sprintf("newGeneticCode(%d,\n\"%s\",\n\"%s\",\n\"%s\",\n\"%s\")",
		gc.Id, gc.Name, gc.ShortName, gc.Ncbieaa, gc.Sncbieaa)
}

func Unqoute(s string) (string, error) {
	if (!strings.HasPrefix(s, "\"")) || (!strings.HasSuffix(s, "\"")) {
		return "", errors.New("String is not quoted")
	}
	return strings.Trim(s, "\""), nil
}

func ANumMinus(b byte) bool {
	r := rune(b)
	return r == '-' || unicode.IsLetter(r) || unicode.IsDigit(r) 
}

func Asn1Split(data []byte, atEOF bool) (int, []byte, error) {
	i := 0
	
	for ; i < len(data) ; i ++ {
		if !unicode.IsSpace(rune(data[i])) {
			break
		}
	}
	data = data[i:]
	advance := i

	if len(data) == 0 {
		return 0, nil, nil
	}

	switch data[0] {
	case '-':
		if len(data) < 2 {
			if atEOF {
				return 0, nil, errors.New("unexpected end of file")
			}
			return advance, nil, nil
		}
		if data[1] == '-' {
			a, t, e := bufio.ScanLines(data, atEOF)
			return a + advance, t, e
		} else {
			return 0, nil, errors.New("unexpected character after '-'")
		}
	case ':':
		if len(data) < 3 {
			if atEOF {
				return 0, nil, errors.New("unexpected end of file")
			}
			return advance, nil, nil
		}
		if data[1] == ':' && data[2] == '=' {
			return advance + 3, data[:3], nil
		} else {
			return 0, nil, errors.New("unxepected character after ':'")
		}
	case '"':
		for i := 1; i < len(data); i++ {
			if data[i] == '"' {
				return advance + i + 1, data[:i+1], nil
			}
		}
		if !atEOF {
			return advance, nil, nil
		} else {
			return 0, nil, errors.New("unfinished string literal")
		}
	case '{', '}', ',':
		return advance + 1, data[:1], nil
	}
	if ANumMinus(data[0]) {
		i := 1
		for ; i < len(data); i++ {
			if !ANumMinus(data[i]) {
				break
			}
		}
		if i == len(data) {
			if atEOF {
				return advance + i, data, nil
			}
			return advance, nil, nil
		}
		return advance + i, data[:i], nil
	}
	return 0, nil, errors.New("Unknown token")
}

func ParseAsn1(rd io.Reader) (res []GeneticCode, err error) {
	scanner := bufio.NewScanner(rd)

	scanner.Split(Asn1Split)
	
	mode := NORMAL

	var gc GeneticCode
	var parName string

	for scanner.Scan() {
		text := scanner.Text()
		if strings.HasPrefix(text, "--") {
			continue
		}
		
		switch mode {
		case NORMAL:
			if text != "Genetic-code-table" {
				return nil, errors.New("Expecting 'Genetic-code-table'")
			}
			mode = TABLE
		case TABLE:
			if text != "::=" {
				return nil, errors.New("Expecting '::='")
			}
			mode = ASSIGN
		case ASSIGN:
			if text != "{" {
				return nil, errors.New("Expecting '{'")
			}
			mode = LIST
		case LIST:
			switch text {
			case "{":
				gc = GeneticCode{}
				mode = ELEMENT
			case "}":
				mode = END
			default:
				return nil, errors.New("Expecting '{' or '}'")
			}
		case ELEMENT:
			parName = text
			mode = ELEMENTPAR
		case ELEMENTPAR:
			switch parName {
			case "name":
				uq, err := Unqoute(text)
				if err != nil {
					return nil, err
				}
				uq = strings.Replace(uq, "\n", "", -1)
				if gc.Name == "" {
					gc.Name = uq
				} else {
					gc.ShortName = uq
				}
			case "id":
				id, err := strconv.Atoi(text)
				if err != nil {
					return nil, err
				}
				gc.Id = id
			case "ncbieaa":
				code, err := Unqoute(text)
				if err != nil {
					return nil, err
				}
				gc.Ncbieaa = code
			case "sncbieaa":
				code, err := Unqoute(text)
				if err != nil {
					return nil, err
				}
				gc.Sncbieaa = code
			}
			mode = ELEMENTPRECOMMA
		case ELEMENTPRECOMMA:
			switch text {
			case ",":
				mode = ELEMENT
			case "}":
				res = append(res, gc)
				mode = PRECOMMA
			default:
				return nil, errors.New("Expecting ',' or '}'")
			}
		case PRECOMMA:
			switch text {
			case ",":
				mode = LIST
			case "}":
				mode = END
			default:
				return nil, errors.New("Expecting ',' or '}'")
			}
		case END:
			return nil, errors.New("Unexpected symbols at the end of file")
		}
	}

	err = scanner.Err()
	if err != nil {
		return nil, err
	}

	if mode != END {
		return nil, errors.New("Unexpected end of stream")
	}

	return
}

func main() {
	if len(os.Args) < 2 {
		fmt.Println("please specify a gc file in asn1 format")
		os.Exit(1)
	}
	f, err := os.Open(os.Args[1])
	if err != nil {
		fmt.Println(err)
		os.Exit(2)
	}
	gcodes, err := ParseAsn1(f)
	if err != nil {
		fmt.Println(err)
		os.Exit(3)
	}

	fmt.Println("package bio")

	fmt.Println("// GeneticCodes is a map holding genetic codes.")
	fmt.Println("// This file was generated using gcode program from NCBI genetic codes file.")
	fmt.Println("var GeneticCodes = map[int]*GeneticCode{")
	for _, gc := range gcodes {
		fmt.Printf("%d: %s,\n", gc.Id, gc.GoString())
	}
	fmt.Println("}")
}
