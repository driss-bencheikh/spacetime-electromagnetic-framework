% esframework.bst -- bibliography style of "The Spacetime-Electromagnetic Framework".
% Writes each entry as   Author, \emph{Title}, <publication data>.   in the order given by the field sortkey.
% Article:  journal volume (year) pages note (with journalemph: Author, Title, \emph{Journal} ...);  book:  publisher, year;
% misc:  howpublished.
ENTRY { author title journal volume year pages note publisher howpublished journalsep journalemph tailsep sortkey } {} {}

FUNCTION {head}
{ author empty$
    { "\emph{" title * "}" * }
    { author ", \emph{" * title * "}" * }
  if$
}

FUNCTION {bibitem}
{ newline$ "\bibitem{" cite$ * "}" * write$ newline$ }

FUNCTION {finish}
{ "." * write$ newline$ }

FUNCTION {article}
{ bibitem
  journalemph empty$
    { head ", " * journal * journalsep empty$ { " " * } { "~" * } if$ }
    { author ", " * title * ", \emph{" * journal * "} " * }
  if$
  "\textbf{" * volume * "} (" * year * ")" *
  pages empty$ { } { " " * pages * } if$
  note empty$ { } { note * } if$
  finish
}

FUNCTION {book}
{ bibitem
  head ", " * publisher * ", " * year *
  finish
}

FUNCTION {misc}
{ bibitem
  head
  howpublished empty$
    { }
    { tailsep empty$ { ", " * howpublished * } { howpublished * } if$ }
  if$
  finish
}

FUNCTION {default.type} { misc }

READ

FUNCTION {presort} { sortkey 'sort.key$ := }
ITERATE {presort}
SORT

FUNCTION {begin.bib}
{ "\begin{thebibliography}{999}" write$ newline$ }
EXECUTE {begin.bib}
ITERATE {call.type$}
FUNCTION {end.bib}
{ newline$ "\end{thebibliography}" write$ newline$ }
EXECUTE {end.bib}
