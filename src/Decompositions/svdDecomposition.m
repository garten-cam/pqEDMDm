classdef svdDecomposition < pqDecomposition
	%SVDDECOMPOSITION Implements the svd decomposition for the calculation
	% of the matrices.
	methods
    function u = regression(obj, lhs, rhs)
			[Ud,S,V] = svd(rhs,"econ");
			s = diag(S);
			% Efffective rank r of xeval_d
			r = obj.effective_rank(s, rhs);
			Ur = Ud(:,1:r);
			Sr = diag(s(1:r));
			Vr = V(:,1:r); % i.e., Sr\Ur'*rhs it is an orthogonalization
			Dr = Sr\Ur'*lhs;
			% solve it
			u = Vr*Dr;
		end
	end
	methods (Static)
		function r = effective_rank(s, xeval_d)
			% all the singular values greater than an arbitrary number
			r = sum(s > max(size(xeval_d)*eps*s(1)));
		end
	end
end

