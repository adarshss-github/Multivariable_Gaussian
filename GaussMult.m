function [p] = GaussMult(X,Mu,COVM)

%===========================================================================================================================================================
%***********************************************************************************************************************************************************
%
%                                      --------------------------------------------------
%                                                          PREAMBLE
%                                      --------------------------------------------------
%
%[p] = GaussMult(X,Mu,Sig)
%
%By Adarsh S, Ph.D. Candidate IIT Kanpur
%
%Function Description:
%--------------------
% Evaluates the multidimensional Gaussian probability density function
% (pdf) at X
% 
% 
%Input Arguments:
%---------------
%X: Column vector at which pdf is to be evaluated
%Mu: Column vector of mean of the multidimensional Gaussian pdf
%COVM: Covariance matrix of the multidimensional Gaussian pdf
%
%Output Arguments:
%-----------------
%p: The value of the multidimensional Gaussian pdf
%
%
%Example: GaussMult(zeros(4,1),zeros(4,1),eye(4)) ;
%
%
%                                      -----------------------------------------------
%                            *********|| All rights reserved; Adarsh S; May, 2020 || *********
%                                      -----------------------------------------------
%
%
%***********************************************************************************************************************************************************
%===========================================================================================================================================================

N = length(Mu) ;

p = sqrt(1/(2*pi)^N/det(COVM))*exp( -0.5*(X-Mu)'*(COVM\eye(N))*(X-Mu)) ;


end