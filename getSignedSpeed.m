function [ signedSpeed ] = getSignedSpeed( velocity, orientation )
%calculates the signed velocity based on the scalar product of vectors
%velocity and orientation, which should be row vectors
% velocity - vector, ndim by nsamples
% orientation - vector, ndim by nsamples

if size(velocity,1)>size(velocity,2)
    velocity = velocity';
end

if size(orientation,1)>size(orientation,2)
    orientation = orientation';
end

% normalise orientation
speed = sqrt(sum(velocity.^2));
signedSpeed = sign(sum(velocity.*orientation)).*speed;

end

