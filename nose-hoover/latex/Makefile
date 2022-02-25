# -----------------------------------------------------------------------------
# Set the target macro
ifeq ($(origin target), undefined)
TARGET := main
else ifeq ($(target),)
TARGET := main
else
TARGET := $(target)
endif

SOURCE := $(TARGET).tex
PDFOUT := $(TARGET).pdf


# -----------------------------------------------------------------------------
# Compiler Settings
AR       := ar rcs
RM       := rm -vf
WC       := wc
ECHO     := echo
VIEWER   := open -a Preview
TEXEDIT  := open -a TeXShop
BIBEDIT  := open -a BibDesk
PDFLATEX := /opt/local/bin/pdflatex --interaction batchmode
BIBTEX   := /opt/local/bin/bibtex


# -----------------------------------------------------------------------------
# Rules
.PHONY: all
all: pdf

.PHONY: pdf
pdf: $(PDFOUT)
	$(VIEWER) $(PDFOUT)

.PHONY: edit
edit: $(SOURCE)
	$(TEXEDIT) $(SOURCE)

.PHONY: clean
clean:
	$(RM) $(PDFOUT) *.log *.aux *.toc *.bbl *.blg *.out

$(PDFOUT):
	$(PDFLATEX) $(SOURCE)
	$(BIBTEX) $(TARGET)
	$(PDFLATEX) $(SOURCE)
	$(PDFLATEX) $(SOURCE)
