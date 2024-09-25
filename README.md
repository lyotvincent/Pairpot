# Pairpot: A Database with Real-time Lasso-based Analysis Tailored for Paired Single-cell and Spatial Transcriptomics

Paired single-cell and spatial resolved transcriptomics (SRT) data supplement each other, providing in-depth insights into biological processes and disease mechanisms.
Previous SRT databases have limitations in curating sufficient single-cell and SRT pairs (SC-SP pairs) and providing real-time heuristic analysis, which hinder the effort to uncover potential biological insights.
Here, we developed Pairpot (http://pairpot.bioxai.cn), a database tailored for paired single-cell and SRT data with real-time heuristic analysis. 
Pairpot curates 99 high-quality pairs including 1,425,656 spots from 299 datasets, and creates the association networks.
It constructs the curated pairs by integrating multiple slices and establishing potential associations between single-cell and SRT data. 
On this basis, Pairpot adopts semi-supervised learning that enables  real-time heuristic analysis for SC-SP pairs where Lasso-View refines the user-selected SRT domains within milliseconds, Pair-View infers cell proportions of spots based on user-selected cell types in real-time, and Layer-View displays SRT slices using a 3D hierarchical layout. 
Experiments demonstrated Pairpot's efficiency in identifying heterogeneous domains and cell proportions. 

## FrontEnd

This project was bootstrapped with [Create React App](https://github.com/facebook/create-react-app).


Built using:

- Front-end library: React
- CSS framework: React-bootstrap
- CSS animations library: Animate.css

In the /frontend, you can run:

### `npm start`

Runs the front-end in the development mode.\
Open [http://localhost:6634](http://localhost:6634) to view it in your browser.

The page will reload when you make changes.\
You may also see any lint errors in the console.

### `npm test`

Launches the test runner in the interactive watch mode.\
See the section about [running tests](https://facebook.github.io/create-react-app/docs/running-tests) for more information.

### `npm run build`

Builds the app for production to the `build` folder.\
It correctly bundles React in production mode and optimizes the build for the best performance.

The build is minified and the filenames include the hashes.\
Your app is ready to be deployed!

See the section about [deployment](https://facebook.github.io/create-react-app/docs/deployment) for more information.

## Back-End

This back-ENd was powered by flask and sqlite via python.

In the /backend you can run:

## `python main.py`

Runs the back-end in the development mode.\
Open [http://localhost:5522](http://localhost:5522) to access the api for the browser.

