std::unique_ptr<CModel> CModel::gradientStep(const Tree & tree, const StateMap & obs, ctpl::thread_pool & threadPool,std::mt19937 & gen, sim::SPolicy & simPolicy) const {
    //Establish the gradient
    //Note A regularized gradient with a magnitude of one is returned
    GradientMap gradient = this->estimateGradient(tree,obs,threadPool,gen,simPolicy);
    //for(const auto & pair : gradient){
    //fprintf(stderr,"%s:%0.4f\t",pair.first.c_str(),pair.second);
    //}
    //fprintf(stderr,"\n");
    //Calculate the magnitude of the gradient vector
     
    //Use {a number} of Chebechev spaced points along the gradient line segment
    //  This spacing is the projection onto the line segment of equally spaced points
    //  on a semicircular arc between the ends of the line segment
    //  The ends of the line segment are the current model's parameters and point at
    //  which the gradient line would hit the end of one of the parameter's domain
    //  If all parameters have an infinite domain, we'll use a {multiple} of the
    //  non-zero, smallest parameter's value
    //  If all domains are infinite, and all parameters are at 0, then just go
    //  {nPoints} full gradient steps

    int nPoints = 5; //TODO: Magic Numbers
    int nDegree = 4; //TODO: Magic Numbers
    int maxMultiple = 3; //TODO: Magic Numbers

    //Determine which parameter is the limiter
    double intervalLength = this->getMaxVectorDist(gradient);
    //fprintf(stderr,"intervalLength: %0.04f\n", intervalLength);
    //intervalLength is now the maximum distance to travel in the gradient direction
    
    //Useful paper Barycentric Lagrange Interpolation: Berrut and Trefethen (2004)
    //  10.1137/S0036144502417715
    //Chebyshev points of the second kind on the interval [-1,1] 
    //  x_j = cos((m-j)pi/m)
    //  Rescale onto the inverval [0,intervalLength]:
    //  x_j = (cos(m-j)pi/m) + 1) * intervalLength/2    

    //Evaluate the model at each point
    std::vector<std::pair<double,double>> vPoints = {std::make_pair(0.0,this->getNLogP())};
    for(int j = 1; j < nPoints - 1; j++){
        ParamMap newParam = this->parameters;
        double stepSize = (std::cos((nPoints - 1 - j) * stats::pi / double(nPoints - 1)) + 1)/2;
        //fprintf(stderr, "%d, stepSize%0.04f\n",j,stepSize);
        stepSize *= intervalLength;
        for(auto & pair : newParam){
            //The stepSize for each parameter is the cosine of the angle between 
            //the gradient and the parameter's basis vector
            //cos(θ_i) = slope_i / gradientMagnitude(1)
            //A single step along the gradient translates to a cos(θ_i) step along the
            //bais :: slope_i / |Δ|
            //We want to minimize the function so we move in the opposite direction of
            //the gradient (-=)
            pair.second.value -= stepSize * gradient.at(pair.first);
        } 
        //fprintf(stderr,"ConstructModel\n");
        std::unique_ptr<CModel> pointModel = this->constructAdjacentModel(newParam,0);
        //std::cerr << *pointModel << "\n";
        pointModel->evaluate(tree,obs,threadPool,gen,simPolicy);
        //fprintf(stderr,"EvaluatedModel\n");
        vPoints.push_back(std::make_pair(stepSize,pointModel->getNLogP()));
    }

    //Infer the approximate location of a minimum as between points where the sign
    //changed from negative to positive
    //If there are multiple go with the pair of points with the lowest total value
    //If there is no flipping of signs go to the last pair
    //Pairs are indexed by the index of the first point
    auto bestPointPair = std::make_pair(nPoints-2,vPoints[nPoints-2].second+vPoints[nPoints-1].second);
    //Determine the sign of change between points
    int lastSign = -1;
    for(int j = 0; j < vPoints.size() -1; j++){
        double diff = vPoints[j+1].second - vPoints[j].second;
        double total = 2.0*vPoints[j].second + diff;
        int sign = (diff < 0.0)  ? -1 : 1;
        if(sign > 0 && lastSign < 0){ //Found a flip
            if(j < bestPointPair.first || total < bestPointPair.second){
                bestPointPair = std::make_pair(j,total);
            }
        }
        lastSign = sign;
    }
    //bestPointPair.first is now the index of a point which is lower than the points on
    //either side (Or is the second to last point; Could be the first point which
    //doesn't have a point before it)
    //Set the final Step size to be average of itself and the previous point weighted
    //by the 'slopes' on either side
    double rightSlopeEst = (vPoints[bestPointPair.first+1].second - vPoints[bestPointPair.first].second) / (vPoints[bestPointPair.first+1].first - vPoints[bestPointPair.first].first);
    double leftSlopeEst = rightSlopeEst; //Default in case of first point index
    if(bestPointPair.first > 0){
        (vPoints[bestPointPair.first].second - vPoints[bestPointPair.first-1].second) / (vPoints[bestPointPair.first].first - vPoints[bestPointPair.first-1].first);
    }
    double totalSlope = std::abs(rightSlopeEst) + std::abs(leftSlopeEst);
    double bestStep = std::abs(leftSlopeEst) * vPoints[bestPointPair.first].first / totalSlope
                + std::abs(rightSlopeEst) * vPoints[bestPointPair.first+1].first / totalSlope;
    //fprintf(stderr,"BestStep: %0.04f\n",bestStep);

    //Construct the parameters at the best Step
    ParamMap newParam = this->parameters;
    for(auto & pair : newParam){
        pair.second.value -= bestStep * gradient.at(pair.first);
    }
    //We are going to say that the probability of proposing this new model is the same
    //as the parent model, as that choice was random, but from there it is
    //deterministic
    return this->constructAdjacentModel(newParam,this->logHastingsRatio);
}
