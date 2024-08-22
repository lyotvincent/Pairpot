// import { useState, useEffect } from 'react'
import contactImg from "../assets/img/contact.png"
import {ConfigProvider, Divider, Layout, Space } from 'antd'
import 'animate.css'
import { Card, Row, Col } from 'antd'
// import { useLocation, useNavigate } from 'react-router-dom'
// import Loading from './charts/Loading'

// steps:
import browser1Img from "../assets/img/steps/browser1.png"
import browser2Img from "../assets/img/steps/browser2.png"
import browser3Img from "../assets/img/steps/browser3.png"
import browser4Img from "../assets/img/steps/browser4.png"

import pair1Img from "../assets/img/steps/pair1.png"
import pair2Img from "../assets/img/steps/pair2.png"
import pair3Img from "../assets/img/steps/pair3.png"

import lasso1Img from "../assets/img/steps/lasso1.png"
import lasso2Img from "../assets/img/steps/lasso2.png"
import lasso3Img from "../assets/img/steps/lasso3.png"
import lasso4Img from "../assets/img/steps/lasso4.png"
import lasso5Img from "../assets/img/steps/lasso5.png"
import lasso6Img from "../assets/img/steps/lasso6.png"





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
              Data Search and Browse
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
                style={{ width: '50%', margin: '2%',marginLeft: '2%', marginRight: '1%', paddingLeft: '1rem', paddingRight: '1rem' }}
                title={<h3 style={{ marginTop: "1rem" }}> Usage1</h3>}
                cover={<img src={browser1Img} alt="Browse-step1"/>}>
                <Card.Meta classNames='header' description={"Users can search for datasets through any words of interest, or apply filters through preset keywords in species, tissues, technologies, and topics."}></Card.Meta>
              </Card>
              <Card hoverable
                style={{ width: '44%', margin: '2%',marginLeft: '1%', marginRight: '2%', paddingLeft: '1rem', paddingRight: '1rem' }}
                title={<h3 style={{ marginTop: "1rem" }}> Usage2</h3>}
                cover={<img src={browser2Img} alt="Browse-step2"/>}>
                <Card.Meta classNames='header' description={"In the extensive list, each item contains a title, summary, pair scores, and action buttons. By clicking ’Visualization’, users can access the analysis results, and perform lasso-based real-time analysis. By clicking ’Descriptions’, users can access detailed meta-information."}></Card.Meta>
              </Card>
              </Row>
              <Row>
              <Card hoverable
                style={{ width: '50%', margin: '2%',marginLeft: '2%', marginRight: '1%', paddingLeft: '1rem', paddingRight: '1rem' }}
                title={<h3 style={{ marginTop: "1rem" }}> Usage3</h3>}
                cover={<img src={browser3Img} alt="Browse-step3"/>}>
                <Card.Meta classNames='header' description={"In the extensive table, Users can browse more items. Actions are the same as the extensive list."}></Card.Meta>
              </Card>
              <Card hoverable
                style={{ width: '44%', margin: '2%',marginLeft: '1%', marginRight: '2%', paddingLeft: '1rem', paddingRight: '1rem' }}
                title={<h3 style={{ marginTop: "1rem" }}> Usage4</h3>}
                cover={<img src={browser4Img} alt="Browse-step4"/>}>
                <Card.Meta classNames='header' description={"Each node denotes a dataset and its color denotes the dataset’s species. Click ‘Display All’ to display all SC-SP pairs. Click ‘Display SRT Rings’ display nodes with self-loops in the network. Click ‘Display SC-SP pairs’ to display no-loops networks."}></Card.Meta>
              </Card>
              </Row>
              </Col>
            <Divider/>
            <Row>
            <Col span={20} offset={1}>
            <div className='headline'>
              <h2>
                {"Lasso View"}
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
                style={{ width: '36%', margin: '2%',marginLeft: '2%', marginRight: '1%', paddingLeft: '1rem', paddingRight: '1rem' }}
                title={<h3 style={{ marginTop: "1rem" }}> Step1</h3>}
                cover={<img src={lasso1Img} alt="Lasso-View-step1"/>}>
                <Card.Meta classNames='header' description={"Click the Lasso Tools at TopRight, then select your cells of interests in the left chart. You can change the select mode through `Settings` button. "}></Card.Meta>
              </Card>
              <Card hoverable
                style={{ width: '28%', margin: '2%', marginLeft: '1%', marginRight: '1%', paddingLeft: '1rem', paddingRight: '1rem' }}
                title={<h3 style={{marginTop: "1rem"}}>Step2</h3>}
                cover={<img src={lasso2Img} alt="Lasso-View-step2"></img>}>
                <Card.Meta classNames='header' description={"Click the `Rename` button and type in the name of your selected cells. Your annotations are in the legends (bottom of the charts)."}></Card.Meta>
              </Card>
              <Card hoverable
                style={{ width: '28%', margin: '2%', marginLeft: '1%', marginRight: '2%', paddingLeft: '1rem', paddingRight: '1rem' }}
                title={<h3 style={{ marginTop: "1rem"}}>Step3</h3>}
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
                style={{ width: '40%', margin: '2%',marginLeft: '1%',marginRight: '1%', paddingLeft: '1rem', paddingRight: '1rem' }}
                title={<h3 style={{ marginTop: "1rem", padding: 0 }}> Step5</h3>}
                cover={<img src={lasso5Img} alt="Lasso-View-step5"></img>}>
                <Card.Meta classNames='header' description={"Lasso your new anntations as same as step 1. You can delete your annotations through `Delete` button."}></Card.Meta>
              </Card>

              <Card hoverable
                style={{ width: '26%', margin: '2%',marginLeft: '1%', paddingLeft: '1rem', paddingRight: '1rem' }}
                title={<h3 style={{ marginTop: "1rem", padding: 0 }}>Step6</h3>}
                cover={<img src={lasso6Img} alt="Lasso-View-step6"></img>}>
                <Card.Meta classNames='header' description={"Save your annotations as a JSON file."}></Card.Meta>
              </Card>
            </Row>
            </Col>

            <Divider/>
            <Row>
              <Col span={20} offset={1}>
            <div className='headline'>
              <h2>
                {"Pair View"}
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
