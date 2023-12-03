#define DEBUG_MODULE "deckPinSwitch"
#include "debug.h"

#include "deck.h"

#include "param.h"
#include "stm32fxxx.h"
#include "FreeRTOS.h"
#include "timers.h"

static bool isInit;

#define TIM_PERIF RCC_APB1Periph_TIM5
#define TIM TIM5
#define TIM_DBG DBGMCU_TIM5_STOP
#define TIM_SETCOMPARE TIM_SetCompare2
#define TIM_GETCAPTURE TIM_GetCapture2

#define GPIO_POS_PERIF RCC_AHB1Periph_GPIOA
#define GPIO_POS_PORT GPIOA
#define GPIO_POS_PIN GPIO_Pin_2 // TIM5_CH3
#define GPIO_AF_POS_PIN GPIO_PinSource2
#define GPIO_AF_POS GPIO_AF_TIM5

#define PERIOD 0x40

static void setUpHwTimer()
{
    // Init structures
    GPIO_InitTypeDef GPIO_InitStructure;
    TIM_TimeBaseInitTypeDef TIM_TimeBaseStructure;
    TIM_OCInitTypeDef TIM_OCInitStructure;

    // Clock the gpio and the timers
    RCC_AHB1PeriphClockCmd(GPIO_POS_PERIF, ENABLE);
    RCC_APB1PeriphClockCmd(TIM_PERIF, ENABLE);

    // Configure the GPIO for the timer output
    GPIO_InitStructure.GPIO_Mode = GPIO_Mode_AF;
    GPIO_InitStructure.GPIO_OType = GPIO_OType_PP;
    GPIO_InitStructure.GPIO_PuPd = GPIO_PuPd_DOWN;
    GPIO_InitStructure.GPIO_Speed = GPIO_Speed_2MHz;
    GPIO_InitStructure.GPIO_Pin = GPIO_POS_PIN;
    GPIO_Init(GPIO_POS_PORT, &GPIO_InitStructure);

    // Map timers to alternate functions
    GPIO_PinAFConfig(GPIO_POS_PORT, GPIO_AF_POS_PIN, GPIO_AF_POS);

    // Timer configuration
    TIM_TimeBaseStructure.TIM_Period = PERIOD - 1;
    TIM_TimeBaseStructure.TIM_Prescaler = 0x05;
    TIM_TimeBaseStructure.TIM_ClockDivision = TIM_CKD_DIV2;
    TIM_TimeBaseStructure.TIM_CounterMode = TIM_CounterMode_Up;
    TIM_TimeBaseStructure.TIM_RepetitionCounter = 0;
    TIM_TimeBaseInit(TIM, &TIM_TimeBaseStructure);

    // PWM channels configuration
    TIM_OCInitStructure.TIM_OCMode = TIM_OCMode_PWM1;
    TIM_OCInitStructure.TIM_OutputState = TIM_OutputState_Enable;
    TIM_OCInitStructure.TIM_Pulse = 0;
    TIM_OCInitStructure.TIM_OCPolarity = TIM_OCPolarity_High;
    TIM_OCInitStructure.TIM_OCIdleState = TIM_OCIdleState_Set;

    // Configure OC3
    TIM_OC3Init(TIM, &TIM_OCInitStructure);
    TIM_OC3PreloadConfig(TIM, TIM_OCPreload_Enable);

    // Enable the timer PWM outputs
    TIM_CtrlPWMOutputs(TIM, ENABLE);
    TIM_SetCompare3(TIM, 0x00);
    TIM_SetCompare4(TIM, 0x00);

    // Enable the timer
    TIM_Cmd(TIM, ENABLE);
}

static void deckPinSwitchInit()
{
    if (isInit)
    {
        return;
    }

    DEBUG_PRINT("*********** Hello Crazyflie 2.1 deck world! from deckPinSwitch*********** \n");

    // HW timer for modulating the LED with > 1MHz
    setUpHwTimer();
    TIM_SetCompare3(TIM, PERIOD / 2);

    isInit = true;
}

static bool deckPinSwitchTest()
{
    DEBUG_PRINT("Hello test passed!from deckPinSwitch\n");
    return true;
}

static const DeckDriver deckPinSwitchDriver = {
    .name = "mydeckPinSwitch",
    .init = deckPinSwitchInit,
    .test = deckPinSwitchTest,

    .usedPeriph = DECK_USING_TIMER5,
    .usedGpio = DECK_USING_PA2,
};

DECK_DRIVER(deckPinSwitchDriver);
