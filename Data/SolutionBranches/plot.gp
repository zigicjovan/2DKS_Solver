set terminal pdfcairo enhanced color font "Helvetica,12"

set output "powerlawK.pdf"
load "plotK.txt"
unset output

set output "powerlawL.pdf"
load "plotL.txt"
unset output

set output "powerlawTK.pdf"
load "plotTK.txt"
unset output

set output "powerlawTL.pdf"
load "plotTL.txt"
unset output

set output "branch.pdf"
load "plotbranch.txt"
unset output
