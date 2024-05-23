import { useAccordionButton } from 'react-bootstrap/AccordionButton'
import PropTypes from 'prop-types'
import { Collapse, ConfigProvider, Space } from 'antd'
import { QuestionCircleOutlined } from '@ant-design/icons'

function CustomToggle({ children, eventKey }) {
  const decoratedOnClick = useAccordionButton(eventKey, () =>
    console.log('totally custom!')
  )

  return (
    <QuestionCircleOutlined
      style={{ border: 'secondary', size: 'xxl' }}
      onClick={decoratedOnClick}></QuestionCircleOutlined>
  )
}

const ToggleAccordion = ({ header, children }) => {
  const items = [
    {
      key: '0',
      children: children,
      label: (
        <Space direction="horizontal" >
          {header}
          <QuestionCircleOutlined
            style={{ size: 'large' }} />
        </Space>
      ),
      style: {
        marginBottom: 0,
      }
    }]
  return (
    <ConfigProvider theme={{
      components:{
        Collapse:{
          headerPadding:0,
          contentPadding:4,
        }
      }
    }}>
      <Collapse ghost items={items} expandIconPosition='end' />
    </ConfigProvider>
  )
}

ToggleAccordion.defaultProps = {
  border: null,
}

export default ToggleAccordion
