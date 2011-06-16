TARGET = ../../tmp/AN-11-148S_temp.pdf
all:	$(TARGET)

$(TARGET): $(wildcard *.tex) $(wildcard *.pdf)
	@cd ../../; ./tdr --style=pas b AN-11-148S
clean:
	@cd ../../; ./tdr clean

final: $(wildcard *.tex) $(wildcard figures/*.pdf)
	@cd ../../; ./tdr --nodraft --style=pas b AN-11-148S
prl: $(wildcard *.tex) $(wildcard figures/*.pdf)
	@cd ../../; ./tdr --nodraft --style=paper,aps,twocolumn b AN-11-148S
