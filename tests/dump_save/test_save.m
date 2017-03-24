(* ::Package:: *)

(*
% $Rev: 1 $
% $Date: 2015-01-29 15:42:05 -0800 (Thu, 29 Jan 2015) $
% $LastChangedBy: jlperla $
% Author: Jesse Perla (c) 2012
% Use, modification and distribution are subject to the 
% Boost Software License, Version 1.0. (See accompanying file 
% LICENSE_ 1_ 0.txt or copy at http://www.boost.org/LICENSE_1_ 0.txt)
*)

notebookDirectoryExists = Check[NotebookDirectory[],0];
outputDirectory = If[notebookDirectoryExists == 0, Directory[], notebookDirectoryExists];
SetDirectory[outputDirectory <> "/log"];

y = x + 2;
y$val = 2.2;
y`scoped = 2.5;

DumpSave["test_save.mx", "Global`"]
