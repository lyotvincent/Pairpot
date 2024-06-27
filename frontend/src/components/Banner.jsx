import { useState, useEffect } from 'react'
import headerImg from '../assets/img/header-img.svg'
import propImg from "../assets/img/database-logo.jpg"
import spImg from "../assets/img/spatial3.png"
import helpImg from "../assets/img/help.png"
import funcImg from "../assets/img/func.png"
import contactImg from "../assets/img/contact.png"
import displayImg from "../assets/img/display.png"
import layerImg from "../assets/img/LayerView.png"
import pairImg from "../assets/img/pair.png"
import lassoImg from "../assets/img/lasso.png"
import netImg from "../assets/img/network.png"
import heatImg from "../assets/img/heatmap.png"
import coverImg from "../assets/img/cover-figure.svg"

import Species from './Species'
import { Layout } from 'antd'
import { ArrowRightCircle } from 'react-bootstrap-icons'
import 'animate.css'
import TrackVisibility from 'react-on-screen'
import { Avatar, Card, Row, Col } from 'antd'
import { SettingOutlined, EditOutlined, EllipsisOutlined } from '@ant-design/icons'
import Landscape3D2 from "./charts/Landscape3D2-1"
const { Meta } = Card

const { Content, Footer } = Layout

const Banner = () => {
  const [loopNum, setLoopNum] = useState(0)
  const [isDeleting, setIsDeleting] = useState(false)
  const [text, setText] = useState('')
  const [delta, setDelta] = useState(300 - Math.random() * 100)
  const [index, setIndex] = useState(1)
  const toRotate = ['Web Developer', 'Web Designer', 'UI/UX Designer']
  const period = 2000
  const contentStyle = {
    height: '160px',
    color: '#fff',
    lineHeight: '50px',
    textAlign: 'center',
    background: '#364d79',
  }


  useEffect(() => {
    let ticker = setInterval(() => {
      tick()
    }, delta)

    return () => {
      clearInterval(ticker)
    }
  }, [text])

  const tick = () => {
    let i = loopNum % toRotate.length
    let fullText = toRotate[i]
    let updatedText = isDeleting
      ? fullText.substring(0, text.length - 1)
      : fullText.substring(0, text.length + 1)

    setText(updatedText)

    if (isDeleting) {
      setDelta((prevDelta) => prevDelta / 2)
    }

    if (!isDeleting && updatedText === fullText) {
      setIsDeleting(true)
      setIndex((prevIndex) => prevIndex - 1)
      setDelta(period)
    } else if (isDeleting && updatedText === '') {
      setIsDeleting(false)
      setLoopNum(loopNum + 1)
      setIndex(1)
      setDelta(500)
    } else {
      setIndex((prevIndex) => prevIndex + 1)
    }
  }

  return (
    <Layout>
      <Content>
        <section className="banner" id="home">
          <Row>
            <Col span={16}>
              <TrackVisibility>
                {({ isVisible }) => (
                  <div
                    className={
                      isVisible ? 'animate__animated animate__fadeIn' : ''
                    }>
                    <div className='headline'>
                      <span className="tagline"></span>
                      <h1>
                        {`Welcome to Pairpot!`}
                        {/* <span
                          className="txt-rotate"
                          dataPeriod="1000"
                          data-rotate='[ "Web Developer", "Web Designer", "UI/UX Designer" ]'>
                          <span className="wrap">{text}</span>
                        </span> */}
                      </h1>
                      <h3>
                        {`A Database and Lasso-Based Platform Tailored for Paired Single-cell and Spatial Transcriptomics`}
                      </h3>
                      <p>
                      Spatial Transcriptomics analysis requires paired single-cell references for cell-type annotation, cell-proportion inference, and cell co-localization detection.
                      However, lacking of paired single-cell and spatial transcriptomics data hinders the in-depth deciphering of tissue architectures and functions. PairPOT curated currently available paired single-cell and spatially resolved transcriptomics data, and performed in-depth analysis including batch effects correction, spatial clustering, markers detection, cell-proportion inference, and cell-cell interaction.
                      </p>
                    </div>
                  </div>
                )}
              </TrackVisibility>
            </Col>
            <Col span={8}>
              <TrackVisibility>
                {({ isVisible }) => (
                  <div
                    className={
                      isVisible ? 'animate__animated animate__zoomIn' : ''
                    }>
                    <img src={coverImg} alt="Header Img" />
                  </div>
                )}
              </TrackVisibility>
            </Col>
          </Row>
          <br />

          <Row gutter={16}>
            {/* 
                      <Row>
            <h2>Web Function</h2>
            <img src={funcImg} style={{width:'4%'}}></img>
          </Row>
            <Col span={8}>
              <Card
                cover={
                  <img
                    alt="example"
                    src={propImg}
                    width="40%"
                  />

                }
                actions={[
                  <SettingOutlined key="setting" />,
                  <EditOutlined key="edit" />,
                  <EllipsisOutlined key="ellipsis" />,
                ]}
              >
                <Meta
                  avatar={<Avatar src="https://api.dicebear.com/7.x/miniavs/svg?seed=8" />}
                  title="Card title"
                  description="This is the description"
                />
              </Card>
            </Col> */}
            {/* <Col span={6}>
              <Card 
              cover={
                <img
                    alt="database"
                    src={propImg}
                    style={{width:'53%',display:'block',margin:'auto'}}
                  // src="https://gw.alipayobjects.com/zos/rmsportal/JiqGstEfoWAOHiTxclqi.png"
                  />
               }
              title="Database" bordered={false}>
                The datasets of single-cell and spatial transcriptomes.
              </Card>
            </Col>
            <Col span={6}>
              <Card 
               cover={
                <img
                    alt="browser"
                    src={spImg}
                    style={{width:'65%',display:'block',margin:'auto'}}
                  // src="https://gw.alipayobjects.com/zos/rmsportal/JiqGstEfoWAOHiTxclqi.png"
                  />
               }
              
              title="Browser" bordered={false}>
                Some cell data correspondence plots.
              </Card>
            </Col> */}
            {/* <Col span={8}>
              <Card 
              cover={
                <img
                    alt="help"
                    src={helpImg}
                    style={{width:'56%',display:'block',margin:'auto'}}
                  // src="https://gw.alipayobjects.com/zos/rmsportal/JiqGstEfoWAOHiTxclqi.png"
                  />
               }
              title="Help" bordered={false}>
                Get in touch with us for question help.
              </Card>
            </Col> */}
          </Row>

          {/* <br/>
          <Row>
          <h2>Some data</h2>
          <img src={displayImg} style={{width:'4%'}}></img>
          </Row>
          <div className='datashow'>
            <Row>
              <div className='datalist'>
                <Col span={25}>
                <Card title="Browser" bordered={false}>
                  <table>
                    <tr>
                      <th>hhhhh</th>
                      <th>Species</th>
                    </tr>
                    <tr>
                      <td>A spatial transcriptome map of developing maize ear</td>
                      <td>Maize</td>
                    </tr>
                    <tr>
                      <td>
                      Spatial transcriptomics map of the embryonic mouse brain – a tool to explore neurogenesis
                      </td>
                      <td>Mus musculus</td>
                    </tr>
                    <tr>
                      <td>SAW: An efficient and accurate data analysis workflow for Stereo-seq spatial transcriptomics</td>
                      <td>Mus musculus</td>
                    </tr>
                    <tr>
                      <td>Spatial transcriptomics combined with single-cell RNA-sequencing unravels the complex infla ...</td>
                      <td>Homo sapiens</td>
                    </tr>
                    <tr>
                      <td>Spatial transcriptomics: The effect of consecutive slices data integration on accurate cell ...</td>
                      <td>Homo sapiens</td>
                    </tr>
                    <tr>
                      <td>Probing pathways by which rhynchophylline modifies sleep using spatial transcriptomics (add ...</td>
                      <td>Mus musculus</td>
                    </tr>
                    <tr>
                      <td>Probing pathways by which rhynchophylline modifies sleep using spatial transcriptomics</td>
                      <td>Mus musculus</td>
                    </tr>
                  </table>
                </Card>
              </Col>
              </div>
              <div className='dataprop'></div>
            </Row>
          </div> */}
              
          <br/>
          <Row>
          <h2>What can Pairpot provide?</h2>
          <img src={displayImg} style={{width:'4%'}}></img>
          </Row>

          <div className='fig-display'>
          <Row>
            <div className='fig-block' style={{width:'43%',margin:'3%'}}>
              <a href="browse#LayerView"><img src={layerImg} ></img></a>
              <h3>Layer View</h3>
              <div className='content'>test: this is the introduction for this fig.</div>
            </div>
            <div className='fig-block'style={{width:'43%',margin:'3%'}}>
              <img src={pairImg} ></img>
              <h3>Pair View</h3>
              <div className='content'>test: this is the introduction for this fig.</div>
            </div>
          </Row>

          <Row>
            <div className='fig-block' style={{width:'35%',margin:'2%'}}>
              <img src={lassoImg}></img>
              <h3>Lasso View</h3>
              <div className='content'>test: this is the introduction for this fig.</div>
            </div>
            <div className='fig-block' style={{width:'23.3%',margin:'2%'}}>
              <img src={netImg}></img>
              <h3>Network</h3>
              <div className='content'>test: this is the introduction for this fig.</div>
            </div>
            <div className='fig-block' style={{width:'26.3%',margin:'2%'}}>
              <img src={heatImg} ></img>
              <h3>Heatmap</h3>
              <div className='content'>test: this is the introduction for this fig.</div>
            </div>
          </Row>

          </div>
          
        </section>
      </Content>

      <Footer>
        <Row>
          <h2>Contact us</h2>
          <img src={contactImg} style={{width:'5%'}}></img>
        </Row>
        <Row>
          <button onClick={() => console.log('connect')}>
            Let’s Connect <ArrowRightCircle size={25} />
          </button>
          <Col offset={3}>
            <p>
              *Zhihan Ruan, Centre for Bioinformatics and Intelligent
              Medicine, Nankai University, rrrzhan@nankai.edu.cn
            </p>
            <p>
              *Jian Liu, Centre for Bioinformatics and Intelligent
              Medicine, Nankai University, jianliu@nankai.edu.cn
            </p>
          </Col>
        </Row>
      </Footer>
    </Layout>
  )
}


export default Banner
