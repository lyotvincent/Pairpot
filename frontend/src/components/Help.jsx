// import { useState, useEffect } from 'react'
import contactImg from "../assets/img/contact.png"
import { ConfigProvider, Divider, Layout, Space } from 'antd'
import 'animate.css'
import { Card, Row, Col } from 'antd'
// import { useLocation, useNavigate } from 'react-router-dom'
// import Loading from './charts/Loading'
import HelpTab from "./utils/HelpTab"
import HelpTab1 from "./utils/HelpTab1"
import HelpTab0 from "./utils/HelpTab0"

// steps:
import browser1Img from "../assets/img/steps/browser1.png"
import browser2Img from "../assets/img/steps/browser2.png"
import browser3Img from "../assets/img/steps/browser3.png"
import browser4Img from "../assets/img/steps/browser4.png"
import browser5Img from "../assets/img/steps/browser5.png"
import browser6Img from "../assets/img/steps/browser6-1.png"

import pair1Img from "../assets/img/steps/pair1.png"
import pair2Img from "../assets/img/steps/pair2.png"
import pair3Img from "../assets/img/steps/pair3.png"

import lasso1Img from "../assets/img/steps/lasso1.png"
import lasso2Img from "../assets/img/steps/lasso2.png"
import lasso3Img from "../assets/img/steps/lasso3.png"
import lasso4Img from "../assets/img/steps/lasso4.png"
import lasso5Img from "../assets/img/steps/lasso5.png"
import lasso6Img from "../assets/img/steps/lasso6.png"

import tool11 from "../assets/img/tool1/1.jpg"
import tool12 from "../assets/img/tool1/2.jpg"
import tool13 from "../assets/img/tool1/3.jpg"
import tool14 from "../assets/img/tool1/4.jpg"
import tool15 from "../assets/img/tool1/5.jpg"
import tool16 from "../assets/img/tool1/6.jpg"

import tool21 from "../assets/img/tool2/0.jpg"
import tool22 from "../assets/img/tool2/1.jpg"
import tool23 from "../assets/img/tool2/2.jpg"
import tool24 from "../assets/img/tool2/3.jpg"






// const { enterLoading, quitLoading } = Loading
// const { Meta } = Card

const { Content, Footer } = Layout

const Help = () => {
  // const [loopNum, setLoopNum] = useState(0)
  // const [isDeleting, setIsDeleting] = useState(false)
  // const [text, setText] = useState('')
  // const [delta, setDelta] = useState(300 - Math.random() * 100)
  // const [index, setIndex] = useState(1)
  // const toRotate = ['Web Developer', 'Web Designer', 'UI/UX Designer']
  // const period = 2000
  // const navigate = useNavigate()
  // const [loading, setLoading] = useState([])
  // const [loadText, setLoadText] = useState("Try an example >>>")
  // const contentStyle = {
  //   height: '160px',
  //   color: '#fff',
  //   lineHeight: '50px',
  //   textAlign: 'center',
  //   background: '#364d79',
  // }


  // useEffect(() => {
  //   let ticker = setInterval(() => {
  //     tick()
  //   }, delta)

  //   return () => {
  //     clearInterval(ticker)
  //   }
  // }, [text])

  // const tick = () => {
  //   let i = loopNum % toRotate.length
  //   let fullText = toRotate[i]
  //   let updatedText = isDeleting
  //     ? fullText.substring(0, text.length - 1)
  //     : fullText.substring(0, text.length + 1)

  //   setText(updatedText)

  //   if (isDeleting) {
  //     setDelta((prevDelta) => prevDelta / 2)
  //   }

  //   if (!isDeleting && updatedText === fullText) {
  //     setIsDeleting(true)
  //     setIndex((prevIndex) => prevIndex - 1)
  //     setDelta(period)
  //   } else if (isDeleting && updatedText === '') {
  //     setIsDeleting(false)
  //     setLoopNum(loopNum + 1)
  //     setIndex(1)
  //     setDelta(500)
  //   } else {
  //     setIndex((prevIndex) => prevIndex + 1)
  //   }
  // }

  return (
    <ConfigProvider theme={{
      components: {
        Card: {
          colorBgContainer: "#fff",
          headerFontSize: 32,
          headerHeight: 48,
          borderRadiusLG: 12,
          colorTextDescription: "rgba(0, 0, 0, 0.6)",
        }
      }
    }}>
      <Layout>
        <Content>
          <section className="banner" id="help">
            <Row>
              <Col span={20} offset={1}>
                <div className="headline">
                  <h2>
                    1. Data Search and Browse
                  </h2>
                  <h3>
                    Enables efficient navigation through the collection of datasets.
                  </h3>
                  <p>
                    Users can browse through an extensive list, table, and association
                    networks of all single-cell and SRT datasets. After filtering by keywords, users can take an
                    overview of statistics and word-clouds about the filtered results.
                  </p>
                </div>
              </Col>
            </Row>
            <Col span={22} offset={1}>
              <Row>
                <Card hoverable
                  style={{ width: '50%', margin: '2%', marginLeft: '2%', marginRight: '1%', paddingLeft: '1rem', paddingRight: '1rem' }}
                  title={<h3 style={{ marginTop: "1rem" }}> Usage1</h3>}
                  cover={<img src={browser1Img} alt="Browse-step1" />}>
                  <Card.Meta classNames='header' description={"Users can search for datasets through any words of interest, or apply filters through preset keywords in species, tissues, technologies, and topics."}></Card.Meta>
                </Card>
                <Card hoverable
                  style={{ width: '44%', margin: '2%', marginLeft: '1%', marginRight: '2%', paddingLeft: '1rem', paddingRight: '1rem' }}
                  title={<h3 style={{ marginTop: "1rem" }}> Usage2</h3>}
                  cover={<img src={browser2Img} alt="Browse-step2" />}>
                  <Card.Meta classNames='header' description={"In the extensive list, each item contains a title, summary, pair scores, and action buttons. By clicking ’Visualization’, users can access the analysis results, and perform lasso-based real-time analysis. By clicking ’Descriptions’, users can access detailed meta-information."}></Card.Meta>
                </Card>
              </Row>
              <Row>
                <Card hoverable
                  style={{ width: '50%', margin: '2%', marginLeft: '2%', marginRight: '1%', paddingLeft: '1rem', paddingRight: '1rem' }}
                  title={<h3 style={{ marginTop: "1rem" }}> Usage3</h3>}
                  cover={<img src={browser3Img} alt="Browse-step3" />}>
                  <Card.Meta classNames='header' description={"In the extensive table, Users can browse more items. Actions are the same as the extensive list."}></Card.Meta>
                </Card>
                <Card hoverable
                  style={{ width: '44%', margin: '2%', marginLeft: '1%', marginRight: '2%', paddingLeft: '1rem', paddingRight: '1rem' }}
                  title={<h3 style={{ marginTop: "1rem" }}> Usage4</h3>}
                  cover={<img src={browser4Img} alt="Browse-step4" />}>
                  <Card.Meta classNames='header' description={"Each node denotes a dataset and its color denotes the dataset’s species. Click ‘Display All’ to display all SC-SP pairs. Click ‘Display SRT Rings’ display nodes with self-loops in the network. Click ‘Display SC-SP pairs’ to display no-loops networks."}></Card.Meta>
                </Card>
              </Row>
              <Row>
                <Card hoverable
                  style={{ width: '50%', margin: '2%', marginLeft: '2%', marginRight: '1%', paddingLeft: '1rem', paddingRight: '1rem' }}
                  title={<h3 style={{ marginTop: "1rem" }}> Usage5</h3>}
                  cover={<img src={browser5Img} alt="Browse-step5" />}>
                  <Card.Meta classNames='header' description={"Users can choose dataset in the search engine to analyse."}></Card.Meta>
                </Card>
                <Card hoverable
                  style={{ width: '44%', margin: '2%', marginLeft: '1%', marginRight: '2%', paddingLeft: '1rem', paddingRight: '1rem' }}
                  title={<h3 style={{ marginTop: "1rem" }}> Usage6</h3>}
                  cover={<img src={browser6Img} alt="Browse-step6" />}>
                  <Card.Meta classNames='header' description={"Click the button of 'Try an example' to analyse data in browse page.Click the interactive graphs and plots with hyperlinks in home page to the Browse section"}></Card.Meta>
                </Card>
              </Row>
            </Col>
            <Divider />
            <Row>
              <Col span={20} offset={1}>
                <div className='headline'>
                  <h2>
                    {"2. Lasso View"}
                  </h2>
                  <h3>
                    Refining the customized cells and domains lassoed by users
                  </h3>
                  <p>
                    Users can select their cells (spots) of interest as a cluster using lasso tools, and then LassoView discovers unselected cells similar to those in the user-defined cluster and eliminates the cells mistakenly selected into the cluster. LassoView is mainly based on label propagation and is optimized by C++, which makes it respond in millisecond-level time.  Operation steps are as following:
                  </p>
                </div>
              </Col>
            </Row>
            <Col span={22} offset={1}>
              <Row>
                <Card hoverable
                  style={{ width: '36%', margin: '2%', marginLeft: '2%', marginRight: '1%', paddingLeft: '1rem', paddingRight: '1rem' }}
                  title={<h3 style={{ marginTop: "1rem" }}> Step1</h3>}
                  cover={<img src={lasso1Img} alt="Lasso-View-step1" />}>
                  <Card.Meta classNames='header' description={"Click the Lasso Tools at TopRight, then select your cells of interests in the left chart. You can change the select mode through `Settings` button. "}></Card.Meta>
                </Card>
                <Card hoverable
                  style={{ width: '28%', margin: '2%', marginLeft: '1%', marginRight: '1%', paddingLeft: '1rem', paddingRight: '1rem' }}
                  title={<h3 style={{ marginTop: "1rem" }}>Step2</h3>}
                  cover={<img src={lasso2Img} alt="Lasso-View-step2"></img>}>
                  <Card.Meta classNames='header' description={"Click the `Rename` button and type in the name of your selected cells. Your annotations are in the legends (bottom of the charts)."}></Card.Meta>
                </Card>
                <Card hoverable
                  style={{ width: '28%', margin: '2%', marginLeft: '1%', marginRight: '2%', paddingLeft: '1rem', paddingRight: '1rem' }}
                  title={<h3 style={{ marginTop: "1rem" }}>Step3</h3>}
                  cover={<img src={lasso3Img} alt="Lasso-View-step3"></img>}>
                  <Card.Meta classNames='header' description={"Click the `Refine` button, choose a refiner (currently `LabelPropagation` is available), and then click `Refine` button to refine your annotations."}></Card.Meta>
                </Card>
              </Row>
            </Col>

            <Col span={22} offset={1}>
              <Row>
                <Card hoverable
                  style={{ width: '26%', margin: '2%', marginLeft: '2%', marginRight: '1%', paddingLeft: '1rem', paddingRight: '1rem' }}
                  title={<h3 style={{ marginTop: "1rem", padding: 0 }}>Step4</h3>}
                  cover={<img src={lasso4Img} alt="Lasso-View-step4"></img>}>
                  <Card.Meta classNames='header' description={"Click `Confirm` button to ensure your current annotation. Then you can lasso for other annotations."}></Card.Meta>
                </Card>

                <Card hoverable
                  style={{ width: '40%', margin: '2%', marginLeft: '1%', marginRight: '1%', paddingLeft: '1rem', paddingRight: '1rem' }}
                  title={<h3 style={{ marginTop: "1rem", padding: 0 }}> Step5</h3>}
                  cover={<img src={lasso5Img} alt="Lasso-View-step5"></img>}>
                  <Card.Meta classNames='header' description={"Lasso your new anntations as same as step 1. You can delete your annotations through `Delete` button."}></Card.Meta>
                </Card>

                <Card hoverable
                  style={{ width: '26%', margin: '2%', marginLeft: '1%', paddingLeft: '1rem', paddingRight: '1rem' }}
                  title={<h3 style={{ marginTop: "1rem", padding: 0 }}>Step6</h3>}
                  cover={<img src={lasso6Img} alt="Lasso-View-step6"></img>}>
                  <Card.Meta classNames='header' description={"Save your annotations as a JSON file."}></Card.Meta>
                </Card>
              </Row>
            </Col>

            <Divider />
            <Row>
              <Col span={20} offset={1}>
                <div className='headline'>
                  <h2>
                    {"3. Pair View"}
                  </h2>
                  <h3>{"Inferring cell proportions of spots using the user-selected cells"}</h3>
                  <p>
                    Pairpot pre-analyzed UCell scores of different cell-type signatures from PanglaoDB for single-cell and SRT data. Users can define customized cell types using lasso tools, and then PairView evaluates the mean UCell scores of user-selected cell types. Finally, PairView generates cell proportions of SRT data by performing NNLS using UCell scores from user-defined cell types and spots in SRT data. Operation steps are as following:
                  </p>
                </div>
              </Col>
            </Row>

            <Col span={22} offset={1}>
              <Row>

                <Card hoverable
                  style={{ width: '24%', margin: '2%', marginLeft: '2%', marginRight: '1%', paddingLeft: '1rem', paddingRight: '1rem' }}
                  title={<h3 style={{ marginTop: "1rem", padding: 0 }}> Step1</h3>}
                  cover={<img src={pair1Img} alt="Pair-View-step1"></img>}>
                  <Card.Meta classNames='header' description={"Click the Lasso Tools at TopRight, then select your cells of interests in the left chart. You can change the select mode through `Settings` button."}></Card.Meta>
                </Card>

                <Card hoverable
                  style={{ width: '44%', margin: '2%', marginLeft: '1%', marginRight: '1%', paddingLeft: '1rem', paddingRight: '1rem' }}
                  title={<h3 style={{ marginTop: "1rem", padding: 0 }}>Step2</h3>}
                  cover={<img src={pair2Img} alt="Pair-View-step2" ></img>}>
                  <Card.Meta classNames='header' description={"Click the `Deconv` buttons to infer cell proportions of spots using the user-selected cells. It would cost a few seconds. You can also refine your selected cell types as same as the operations in LassoView."}></Card.Meta>
                </Card>

                <Card hoverable
                  style={{ width: '24%', margin: '2%', marginLeft: '1%', marginRight: '2%', paddingLeft: '1rem', paddingRight: '1rem' }}
                  title={<h3 style={{ marginTop: "1rem", padding: 0 }}> Step3</h3>}
                  cover={<img src={pair3Img} alt="Pair-View-step3" ></img>}>
                  <Card.Meta classNames='header' description={"Users can drag the bar at the top of charts to adjust the thresholds. Users can also explore cell proportions of different batches through `Batches` in `Settings` button."}></Card.Meta>
                </Card>
              </Row>
            </Col>


            <Divider />
            <Row>
              <Col span={20} offset={1}>
                <div className='headline'>
                  <h2>
                    {"4. Alternative Bioinformatics Tools"}
                  </h2>
                  <h3>
                    Selecting alternative bioinformatics tools by users
                  </h3>
                  <p>
                    Users can browse and select these tools on the website.  Operation steps are as following:
                  </p>
                </div>
              </Col>
            </Row>

            <Row>
              <div style={{ marginLeft: '7rem' }}>
                <h4>Selections of deconvolution tools in Pair-View</h4>
              </div>
            </Row>
            <Col span={22} offset={1}>
              <Row>
                <Card hoverable
                  style={{ width: '36%', margin: '2%', marginLeft: '2%', marginRight: '1%', paddingLeft: '1rem', paddingRight: '1rem' }}
                  title={<h3 style={{ marginTop: "1rem" }}> Step1</h3>}
                  cover={<img src={tool11} alt="tool11" />}>
                  <Card.Meta classNames='header' description={"Click the 'spConfigs' button to see the options for Spatial View. "}></Card.Meta>
                </Card>
                <Card hoverable
                  style={{ width: '28%', margin: '2%', marginLeft: '1%', marginRight: '1%', paddingLeft: '1rem', paddingRight: '1rem' }}
                  title={<h3 style={{ marginTop: "1rem" }}>Step2</h3>}
                  cover={<img src={tool12} alt="tool12"></img>}>
                  <Card.Meta classNames='header' description={"There's some alternative tools in the 'Props', including Cell2Location, Card, RCTD, Seurat, SpaTalk..."}></Card.Meta>
                </Card>
                <Card hoverable
                  style={{ width: '28%', margin: '2%', marginLeft: '1%', marginRight: '2%', paddingLeft: '1rem', paddingRight: '1rem' }}
                  title={<h3 style={{ marginTop: "1rem" }}>Step3</h3>}
                  cover={<img src={tool13} alt="tool13"></img>}>
                  <Card.Meta classNames='header' description={"The default deconvolution tool is Cell2Location, users can select other alternative tools."}></Card.Meta>
                </Card>
              </Row>
              <Row>
                <Card hoverable
                  style={{ width: '26%', margin: '2%', marginLeft: '2%', marginRight: '1%', paddingLeft: '1rem', paddingRight: '1rem' }}
                  title={<h3 style={{ marginTop: "1rem", padding: 0 }}>Step4</h3>}
                  cover={<img src={tool14} alt="tool14"></img>}>
                  <Card.Meta classNames='header' description={"Click the 'apply' button to apply for your settings and wait for the change of the Spatial View on the right."}></Card.Meta>
                </Card>

                <Card hoverable
                  style={{ width: '33%', margin: '2%', marginLeft: '1%', marginRight: '1%', paddingLeft: '1rem', paddingRight: '1rem' }}
                  title={<h3 style={{ marginTop: "1rem", padding: 0 }}> Step5</h3>}
                  cover={<img src={tool15} alt="tool15"></img>}>
                  <Card.Meta classNames='header' description={"Select a cell-type on the Single-cell Data chart."}></Card.Meta>
                </Card>

                <Card hoverable
                  style={{ width: '33%', margin: '2%', marginLeft: '1%', paddingLeft: '1rem', paddingRight: '1rem' }}
                  title={<h3 style={{ marginTop: "1rem", padding: 0 }}>Step6</h3>}
                  cover={<img src={tool16} alt="tool16"></img>}>
                  <Card.Meta classNames='header' description={"The corresponding proportions will be shown on the Spatial Omics Data chart."}></Card.Meta>
                </Card>
              </Row>
            </Col>

            <Row>
              <div style={{ marginLeft: '7rem' }}>
                <h4>Selections of cell/domain interaction tools in Heatmap</h4>
              </div>
            </Row>
            <Col span={22} offset={1}>
              <Row>
                <Card hoverable
                  style={{ width: '31%', margin: '2%', marginLeft: '2%', marginRight: '1%', paddingLeft: '1rem', paddingRight: '1rem' }}
                  title={<h3 style={{ marginTop: "1rem", padding: 0 }}>Step1</h3>}
                  cover={<img src={tool21} alt="heatmap1"></img>}>
                  <Card.Meta classNames='header' description={"Click this button to switch spatial or single-cell dataset for displaying L-R pairs heatmap."}></Card.Meta>
                </Card>
                <Card hoverable
                  style={{ width: '30%', margin: '2%', marginLeft: '1%', marginRight: '1%', paddingLeft: '1rem', paddingRight: '1rem' }}
                  title={<h3 style={{ marginTop: "1rem", padding: 0 }}> Step2</h3>}
                  cover={<img src={tool22} alt="heatmap2"></img>}>
                  <Card.Meta classNames='header' description={"Users can select the alternative tools curresponding to current single-cell or SRT data."}></Card.Meta>
                </Card>

                <Card hoverable
                  style={{ width: '30%', margin: '2%', marginLeft: '1%', paddingLeft: '1rem', paddingRight: '1rem' }}
                  title={<h3 style={{ marginTop: "1rem", padding: 0 }}>Step3</h3>}
                  cover={<img src={tool23} alt="heatmap3"></img>}>
                  <Card.Meta classNames='header' description={"Select a tool and wait for the result."}></Card.Meta>
                </Card>
              </Row>
            </Col>


            <Divider style={{ marginTop: 14, marginBottom: 18 }} />
            <Col span={22} offset={1}>
              <div className="headline">
                <h2>5. How to deploy bioinformatics tools offline?</h2>
              </div>
            </Col>
            <Col span={22} offset={1}>
              <Row style={{ marginTop: "1rem" }}>
                <Card style={{ width: '49%' }} title={"Cell/domain Interaction"}>
                  <HelpTab />
                </Card>
                <Card style={{ width: '49%', marginLeft: '2%' }} title={"Deconvolution"}>
                  <HelpTab0 />
                </Card>
              </Row>
            </Col>
            <Divider style={{ marginTop: 14, marginBottom: 18 }} />
            <Col span={22} offset={1}>
              <div className="headline">
                <h2>6. How to convert anndata to R environment?</h2>
              </div>
            </Col>
            <Col span={22} offset={1}>
              <Card style={{ marginTop: "1rem" }}>
                <HelpTab1 />
              </Card>
            </Col>
          </section>
        </Content>

        <Footer>
          <Row>

          </Row>
          <Row>
            <Col span={4}>
              <Space direction='horizontal'>
                <h3>Contact us
                  {' '}
                  <img src={contactImg} style={{ width: '20%', height: '20%' }} alt="contact"></img>
                </h3>
              </Space>
              <a href="https://beian.miit.gov.cn" target='_blank'
                style={{
                  color: 'rgba(0, 0, 0, 0.6)',
                  textDecoration: 'none',
                  fontSize: 13
                }}>
                津ICP备2024022781号-1
              </a>
            </Col>
            <Col offset={1}>
              <p>
                *Zhihan Ruan, Centre for Bioinformatics and Intelligent
                Medicine, Nankai University, rrrzhan@nankai.edu.cn
              </p>
              <p>
                *Jian Liu, State Key Laboratory of Medical Chemical Biology, College of Computer Science, Nankai University, jianliu@nankai.edu.cn
              </p>
            </Col>
          </Row>
        </Footer>
      </Layout>
    </ConfigProvider>
  )
}


export default Help
