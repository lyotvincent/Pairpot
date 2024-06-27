import { Row, Col } from 'antd'
import { MailchimpForm } from './MailchimpForm'
import logo from '../assets/img/logo.svg'
import navIcon1 from '../assets/img/nav-icon1.svg'
import navIcon2 from '../assets/img/nav-icon2.svg'
import navIcon3 from '../assets/img/nav-icon3.svg'

export const Footer = () => {
  return (
    <footer className="footer">
      <Row>
        <Col span={4} offset={4}>
          <img src={logo} alt="Logo" />
        </Col>
        <Col span={6} sm={6} className="text-center text-sm-end">
          <div className="social-icon">
            <a href="#">
              <img src={navIcon1} alt="Icon" />
            </a>
            <a href="#">
              <img src={navIcon2} alt="Icon" />
            </a>
            <a href="#">
              <img src={navIcon3} alt="Icon" />
            </a>
          </div>
        </Col>
      </Row>
      <Row>
        <Col offset={13}>
          <p>
            Copyright © 2024. Centre for Bioinformatics and Intelligent
            Medicine, Nankai University.
          </p>
        </Col>
      </Row>
    </footer>
  )
}

export default Footer
