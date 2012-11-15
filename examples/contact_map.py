#!/usr/bin/env python
# -*- coding: utf-8 -*-

# linijke coding podaje sie by mozna bylo uzywac np polskich znakow w kodzie

import pycabs
from math import sqrt # importujemy funkcje do liczenia pierwiastka

pdb_model = pycabs.parsePDBfile("../playground/2pcy.pdb") # path to the file in the PDB format
# w pdb_model będą atomy Cα w jednowymiarowej liscie, tzn ze jest dlugosci 3xN i np pdb_model[3] to x-owa wspolrzedna drugiego atomu, a pdb_model[4] to y-kowa wspolrzedna 2-go atomu
print pdb_model # test, czy faktycznie

def distance(i,j,model):
	xi = model[3*i]     # na pozycjach w takiej liscie 3*i beda zawsze x-y
	yi = model[3*i + 1] # a na pozycji 3i + 1 beda zawsze y-ki
	zi = model[3*i + 2]
	
	xj = model[3*j]     # i to samo dla j-tego atomu
	yj = model[3*j+1]   
	zj = model[3*j+2]     
	
	xt = xi - xj # odejmujemy, bo wzor na odleglosc kartezjanska to pierwiastek z (xi - xj)^2 + (yi-yj)^2 + (zi-zj)^2 
	yt = yi - yj
	zt = zi - zj
	
	distance = sqrt(xt*xt + yt*yt + zt*zt)
	return distance

# teraz test (zakomentuj to sobie pozniej - komentuje sie znakiem "#"), zobaczymy czy odleglosci miedzy atomami i-1,i to faktycznie 3.8 Angstroma:

for i in range(1,len(pdb_model)/3): # podzielilem przez trzy bo ta lista ma dlugosc 3N, a teraz i bedzie od 1...N (tak naprawde do 1...N-1)
    odleglosc = distance(i-1,i,pdb_model)
    print "%5.2f" %(odleglosc) # mozna to tak wyswietlic, wtedy bedzie zaokraglalo do 2-go miejsca po przecinku. jak dasz samo print odleglosc to zaokragli do ok 15tego miejsca - brzydko wyglada
    # teraz jak to uruchomisz to wyswietli liste, gdzie bedzie mniej wiecej 3.8
    
# no dobra, a teraz chcemy mape kontaktow
# tworzymy pusta macierz 2x2
from numpy import zeros # importujemy funkcje do tworzenia wygodnych macierzy 2D

l = len(pdb_model)/3 # zapiszemy dlugosc lancucha do zmiennej l, bo uzyjemy jej kilka razy, tak ze zeby nie trzeba bylo za kazdym razem tego liczyc

macierz = zeros((l,l)) # tu ma podwojne nawiasy, bo tak naprawde to (l,l) to jest krotka/tupla. tak jakby jedna zmienna. ona definiuje wymiar macierzy. funkcja nazywa sie zeros, bo tworzy macierz wypelniona zerami

cutoff = 6.0 # definiujemy cutoff kontaktu. jesli odleglosc < cutoff, wtedy mamy kontakt

for i in range(l):
	for j in range(l):
		d = distance(i,j,pdb_model)
		if d <= cutoff:
			macierz[i][j] = 1
		#else:  tego nie trzeba, bo macierz od poczatku ma zera

# i teraz mozemy zapisac taka macierz:
pycabs.heat_map(macierz,"Numer reszty", "Numer reszty","Mapa kontaktow",output_file="kontakty.png")

# powinien sie pojawic plik "kontakty.png" z obrazkiem mapy kontaktow
# no i bedziemy robic podobne rzeczy, tylko dla calej trajektorii, a nie tylko dla jednego modelu
