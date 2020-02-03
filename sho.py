from manimlib.imports import *
from myProjects.fta import Grid, ScreenGrid
import numpy as np

class Setup(Scene):
    CONFIG = {
        "frequency" : 0.5,
        "ceiling_radius" : 3*FRAME_X_RADIUS,
        "n_springs" : 1,
        "amplitude" : 0.6,
        "spring_radius" : 0.15,
    }

    def construct(self):
        self.setup_springs()
        self.wait(5)
        title = TextMobject("What is the maths describing this motion?")
        title.to_edge(UP)
        self.play(ShowCreation(title))
        self.wait(2)


    def setup_springs(self):
        def get_spring(alpha, height = 2):
            t_max = 6.5
            r = self.spring_radius
            s = (height - r)/(t_max**2)
            spring = ParametricFunction(
                lambda t : op.add(
                    r*(np.sin(TAU*t)*UP+np.cos(TAU*t)*LEFT),
                    s*((t_max - t)**2)*RIGHT,
                ),
                t_min = 0, t_max = t_max,
                color = WHITE,
                stroke_width = 2,
            )
            spring.alpha = alpha
            return spring
        alphas = np.linspace(0, 1, self.n_springs)
        bezier([0, 1, 0, 1])
        springs = self.springs = VGroup(*list(map(get_spring, alphas)))


        k_tracker = self.k_tracker = VectorizedPoint()
        t_tracker = self.t_tracker = VectorizedPoint()
        always_shift(t_tracker, RIGHT, 1)
        self.t_tracker_walk = t_tracker
        equilibrium_height = springs.get_height()

        boxes = self.boxes = VGroup()
        box_anims = box_anims = []
        def update_springs(springs):
            for spring in springs:
                k = k_tracker.get_center()[0]
                t = t_tracker.get_center()[0]
                f = self.frequency
                x = spring.get_top()[0]
                A = self.amplitude
                d_height = A*np.cos(TAU*f*t - k*x)
                new_spring = get_spring(spring.alpha, 2+d_height)
                box = Square(side_length = 1)
                box.spring = spring
                box_anim = Mobject.add_updater(
                    box, lambda b: b.move_to(b.spring.get_start()+RIGHT*0.3)
                )
                box.set_fill(opacity = 1)
                box_anim.update(0)
                box_anims.append(box_anim)
                boxes.add(box)

                spring.become(new_spring)

        spring_update_anim = Mobject.add_updater(springs, update_springs)
        self.spring_update_anim = spring_update_anim
        spring_update_anim.update(0)
        self.add(self.t_tracker_walk)
        self.add(*box_anims)
        self.add(spring_update_anim)



class ForcesAndEquations(Scene):
    CONFIG = {
        "frequency" : 0.5,
        "ceiling_radius" : 3*FRAME_X_RADIUS,
        "n_springs" : 1,
        "amplitude" : 0.6,
        "spring_radius" : 0.15,
    }


    def setup_springs(self):
        def get_spring(alpha, height = 2):
            t_max = 6.5
            r = self.spring_radius
            s = (height - r)/(t_max**2)
            spring = ParametricFunction(
                lambda t : op.add(
                    r*(np.sin(TAU*t)*UP+np.cos(TAU*t)*LEFT),
                    s*((t_max - t)**2)*RIGHT,
                ),
                t_min = 0, t_max = t_max,
                color = WHITE,
                stroke_width = 2,
            )
            spring.alpha = alpha
            return spring
        alphas = np.linspace(0, 1, self.n_springs)
        bezier([0, 1, 0, 1])
        springs = self.springs = VGroup(*list(map(get_spring, alphas)))


        k_tracker = self.k_tracker = VectorizedPoint()
        t_tracker = self.t_tracker = VectorizedPoint()
        always_shift(t_tracker, RIGHT, 1)
        self.t_tracker_walk = t_tracker
        equilibrium_height = springs.get_height()

        boxes = self.boxes = VGroup()
        box_anims = box_anims = []
        def update_springs(springs):
            for spring in springs:
                k = k_tracker.get_center()[0]
                t = t_tracker.get_center()[0]
                f = self.frequency
                x = spring.get_top()[0]
                A = self.amplitude
                d_height = A*np.cos(TAU*f*t - k*x)
                new_spring = get_spring(spring.alpha, 2+d_height)
                box = Square(side_length = 1)
                box.spring = spring
                box_anim = Mobject.add_updater(
                    box, lambda b: b.move_to(b.spring.get_start()+RIGHT*0.3)
                )
                box.set_fill(opacity = 1)
                box_anim.update(0)
                box_anims.append(box_anim)
                boxes.add(box)

                spring.become(new_spring)

        spring_update_anim = Mobject.add_updater(springs, update_springs)
        self.spring_update_anim = spring_update_anim
        spring_update_anim.update(0)
        self.add(springs)
        self.add(boxes)

    def construct(self):

        self.setup_springs()

        forceOneText = TexMobject("{F}_{extension} = -kx", color = RED)
        forceTwoText = TexMobject("{F}_{friction} = -b\dot {x} ", color = GREEN)
        forceThreeText = TexMobject("F=ma")
        forceFourText = TexMobject("ma =", " -kx")
        forceFourText[1].set_color(RED)
        forceFiveText = TexMobject("ma =", "-kx", " +",  "(-b\dot{x})")
        forceFiveText[1].set_color(RED)
        forceFiveText[3].set_color(GREEN)
        forceSixText = TexMobject("ma =", " -kx", "- b\dot{x}")
        forceSixText[1].set_color(RED)
        forceSixText[2].set_color(GREEN)
        forceSevenText = TexMobject("m\ddot{x} =", " -kx", "-b\dot{x}")
        forceSevenText[1].set_color(RED)
        forceSevenText[2].set_color(GREEN)
        forceEightText = TexMobject("m\ddot{x}", " + ", "b\dot{x}", " +", " kx", " = 0")
        forceEightText[2].set_color(GREEN)
        forceEightText[4].set_color(RED)
        forceOneArrow = Arrow(set_length = 3, color = RED)
        forceTwoArrow = Arrow(set_length = 1, color = GREEN).scale(0.5)


        forceOneArrow.shift(1.9*RIGHT)
        forceOneArrow.rotate(PI)
        forceTwoArrow.shift(3.4*RIGHT+0.5*DOWN)
        forceOneText.to_corner(UL)
        forceTwoText.move_to(forceOneText.get_center()+RIGHT*4)
        forceThreeText.move_to(forceOneText.get_center()+DOWN)
        forceFourText.move_to(forceThreeText.get_center()+DOWN)
        forceFiveText.move_to(forceFourText.get_center()+DOWN)
        forceSixText.move_to(forceFiveText.get_center())
        forceSevenText.move_to(forceSixText.get_center())
        forceEightText.move_to(forceSevenText.get_center()+DOWN)


        self.play(Write(forceOneText), Write(forceOneArrow))
        self.wait(2)
        self.play(Write(forceThreeText))
        self.wait(2)
        self.play(ReplacementTransform(forceOneText.copy(), forceFourText), ReplacementTransform(forceThreeText.copy(), forceFourText))
        self.wait(2)
        self.play(Write(forceTwoText), Write(forceTwoArrow))
        self.wait(2)
        self.play(ReplacementTransform(forceTwoText.copy(), forceFiveText), ReplacementTransform(forceFourText.copy(), forceFiveText))
        self.wait(2)
        self.play(ReplacementTransform(forceFiveText, forceSixText))
        self.wait(2)
        self.play(ReplacementTransform(forceSixText, forceSevenText))
        self.wait(2)
        self.play(ReplacementTransform(forceSevenText.copy(), forceEightText))
        self.play(Indicate(forceEightText))

        self.wait(2)

class LotsOfRearranging(Scene):
    def construct(self):
        sho_ode = TexMobject("\ddot{x}",  "+",  "\\frac{b}{m}", "\dot{x}", "+", "\\frac{k}{m}", "x",  "= 0")
        sho_ode[0].set_color(RED)
        sho_ode[3].set_color(BLUE)
        sho_ode[6].set_color(GREEN)

        realise = TextMobject("Key Insights:").scale(0.8)
        realiseOne = TextMobject("1) If this was a quadratic, we could solve it").scale(0.8)
        realiseTwo = TextMobject("2) The first and second derivative are proportional to the function").scale(0.8)

        sho_ode.to_edge(UP)
        realise.move_to(realiseOne.get_center()+ UP)
        realiseTwo.move_to(realiseOne.get_center() + DOWN)

        self.play(ShowCreation(sho_ode))
        self.wait(2)
        self.play(ShowCreation(realise))
        self.wait(2)
        self.play(ShowCreation(realiseOne))
        self.wait(2)
        self.play(ShowCreation(realiseTwo))
        self.wait(2)
        self.play(FadeOut(realiseOne), FadeOut(realiseTwo), FadeOut(realise))
        self.wait(2)

        function = TexMobject("x(t) = Ae ^{\\alpha t}")
        function[0].set_color(GREEN)
        disclaimer = TextMobject("After some testing and tinkering with different functions...").scale(0.4)

        disclaimer.move_to(function.get_center() + DOWN)
        self.play(ShowCreation(disclaimer))
        self.wait()
        self.play(ShowCreation(function))
        self.wait()
        self.play(FadeOut(disclaimer))

        self.wait()

        derivativeOne = TexMobject("\dot{x}(t) = \\alpha A e^{\\alpha t}")
        derivativeOne[0].set_color(BLUE)
        derivativeTwo = TexMobject("\ddot{x}(t) = \\alpha^2 A e^{\\alpha t}")
        derivativeTwo[0].set_color(RED)

        derivativeOne.move_to(function.get_center() + DOWN)
        derivativeTwo.move_to(derivativeOne.get_center() + DOWN)

        self.play(ShowCreation(derivativeOne))
        self.wait()
        self.play(ShowCreation(derivativeTwo))
        self.wait(2)

        intermediaryOne = TexMobject("\\alpha^2 A e^{\\alpha t}", "+",  "\\frac{b}{m}", "\dot{x}", "+", "\\frac{k}{m}", "x",  "= 0")
        intermediaryOne.to_edge(UP)
        intermediaryOne[0].set_color(RED)
        intermediaryOne[3].set_color(BLUE)
        intermediaryOne[6].set_color(GREEN)

        self.play(ReplacementTransform(derivativeTwo, intermediaryOne), ReplacementTransform(sho_ode, intermediaryOne))
        self.wait()

        intermediaryTwo = TexMobject("\\alpha^2 A e^{\\alpha t}", "+",  "\\frac{b}{m}", "\\alpha A e^{\\alpha t}", "+", "\\frac{k}{m}", "x",  "= 0")
        intermediaryTwo.to_edge(UP)
        intermediaryTwo[0].set_color(RED)
        intermediaryTwo[3].set_color(BLUE)
        intermediaryTwo[6].set_color(GREEN)

        self.play(ReplacementTransform(intermediaryOne, intermediaryTwo), ReplacementTransform(derivativeOne, intermediaryTwo))
        self.wait()

        intermediaryThree = TexMobject("\\alpha^2 A e^{\\alpha t}", "+",  "\\frac{b}{m}", "\\alpha A e^{\\alpha t}", "+", "\\frac{k}{m}", "Ae ^{\\alpha t}",  "= 0")
        intermediaryThree.to_edge(UP)
        intermediaryThree[0].set_color(RED)
        intermediaryThree[3].set_color(BLUE)
        intermediaryThree[6].set_color(GREEN)

        self.play(ReplacementTransform(intermediaryTwo, intermediaryThree), ReplacementTransform(function, intermediaryThree))
        self.wait(2)

        final = TexMobject("Ae^{\\alpha t} (\\alpha^2+\\frac{b}{m} \\alpha+\\frac{k}{m})=0")
        self.play(ReplacementTransform(intermediaryThree.copy(), final))

        self.wait()

        note = TextMobject("If A = 0, there is no movement - which is boring!")
        note.move_to(final.get_center()+DOWN)
        self.play(ShowCreation(note))

        self.wait()

        quadratic = TexMobject("\\therefore \\alpha^2+\\frac{b}{m} \\alpha+\\frac{k}{m}=0")
        quadratic.move_to(note.get_center()+DOWN)
        self.play(ShowCreation(quadratic))
        self.play(Indicate(quadratic))

        self.wait()

class FindingAlpha(Scene):
    def construct(self):
        quadratic = TexMobject("\\alpha^2 + \\frac{b}{m} \\alpha + \\frac{k}{m}=0").scale(0.7)
        note = TextMobject("Using quadratic formula...").scale(0.7)
        alphaSolution = TexMobject("\\alpha =\\frac{ \\frac{-b}{m} \\pm \\sqrt{ \\frac{b^2}{m^2} - \\frac{4k}{m}} }{2(1)} ").scale(0.7)
        alphaSolutionSimplified = TexMobject("\\alpha =\\frac{-b}{2m} \\pm \\frac{\\sqrt{ \\frac{b^2}{m^2} - \\frac{4k}{m}}}{2}").scale(0.7)
        alphaSolutionFactorised = TexMobject("\\alpha =\\frac { -b }{ 2m } \\pm \\frac { \\sqrt { \\frac { -4k }{ m }  } \\sqrt { 1-\\frac { m }{ 4k } \\frac { b^{ 2 } }{ m^{ 2 } }  }  }{ 2 } ").scale(0.7)
        alphaSolutionFactorisedSimplified = TexMobject("\\alpha =\\frac { -b }{ 2m } \\pm \\frac { 2i \\sqrt{\\frac{k}{m}}\\sqrt { 1- \\frac { b^{ 2 } }{ 4km }  }  }{ 2 } ").scale(0.7)
        omegaPrime = TexMobject("\\omega \' = \\sqrt{\\frac { k }{ m }}  \\sqrt { 1- \\frac { b^{ 2 } }{ 4km }  }").scale(0.5)
        omegaPrime[0].set_color(BLUE)
        finalAlpha = TexMobject("\\alpha = \\frac{-b}{2m} \\pm i", "\\omega\'").scale(0.7)
        finalAlpha[1].set_color(BLUE)

        quadratic.to_edge(UP)
        note.move_to(quadratic.get_center()+DOWN*1.3)
        alphaSolution.move_to(quadratic.get_center()+DOWN*1.3)
        alphaSolutionSimplified.move_to(alphaSolution.get_center()+DOWN*1.3)
        alphaSolutionFactorised.move_to(alphaSolutionSimplified.get_center()+DOWN*1.3)
        alphaSolutionFactorisedSimplified.move_to(alphaSolutionFactorised.get_center()+DOWN*1.3)
        omegaPrime.move_to(quadratic.get_center()+RIGHT*3)
        finalAlpha.move_to(alphaSolutionFactorisedSimplified.get_center()+DOWN*1.3)

        self.play(Write(quadratic))
        self.wait(2)
        self.play(Write(note))
        self.play(FadeOut(note))
        self.wait(2)
        self.play(ReplacementTransform(quadratic.copy(), alphaSolution))
        self.wait(2)
        self.play(ReplacementTransform(alphaSolution.copy(), alphaSolutionSimplified))
        self.wait(2)
        self.play(ReplacementTransform(alphaSolutionSimplified.copy(), alphaSolutionFactorised))
        self.wait(2)
        self.play(Indicate(alphaSolutionFactorised, color = RED))
        self.wait(2)
        self.play(ReplacementTransform(alphaSolutionFactorised.copy(), alphaSolutionFactorisedSimplified))
        self.wait(2)
        self.play(Write(omegaPrime))
        self.wait(2)
        self.play(ReplacementTransform(alphaSolutionFactorisedSimplified.copy(), finalAlpha), ReplacementTransform(omegaPrime.copy(), finalAlpha))
        self.wait(2)
        self.play(Indicate(finalAlpha))
        self.wait()


class FindingDisplacement(Scene):
    def construct(self):
        function = TexMobject("x(t) = ","A", "e ^ ", "{\\alpha", " t}").scale(0.7)
        function[3].set_color(RED)
        alpha = TexMobject("\\alpha = \\frac{-b}{2m} + ", "i\\omega\'").scale(0.7)
        alpha[0].set_color(RED)
        alpha[1].set_color(BLUE)
        omegaPrime = TexMobject("\\omega \' = \\sqrt{\\frac { k }{ m }}  \\sqrt { 1- \\frac { b^{ 2 } }{ 4km }  }").scale(0.7)
        omegaPrime.set_color(BLUE)

        function.to_corner(UL)
        alpha.to_edge(UP)
        omegaPrime.to_corner(UR)

        self.play(Write(function), Write(alpha), Write(omegaPrime))
        self.wait(2)
        combinedFunction = TexMobject("x(t) = Ae^", "{(\\frac{-b}{2m}+", "i\\omega\'", " )t}")
        combinedFunction[1].set_color(RED)
        combinedFunction[2].set_color(BLUE)
        splitFunction = TexMobject("x(t) = ", "A", "e^{\\frac{-b}{2m}t}", "e^{i\\omega\'t}")
        splitFunction[2].set_color(RED)
        splitFunction[3].set_color(BLUE)

        self.play(ReplacementTransform(function.copy(), combinedFunction), ReplacementTransform(alpha.copy(), combinedFunction))
        self.wait(2)
        self.play(ReplacementTransform(combinedFunction, splitFunction))
        self.wait(2)

        functionUsed = TexMobject("x(t) = ", "A", "e^{\\frac{-b}{2m}t}", "e^{i\\omega\'t}")
        functionUsed[2].set_color(RED)
        functionUsed[3].set_color(BLUE)
        functionUsed.to_edge(UP)
        self.play(FadeOut(function), FadeOut(alpha), FadeOut(omegaPrime), ReplacementTransform(splitFunction, functionUsed))


        A = 2
        b = 0.4
        m = 1
        k = 7
        phi = 0
        ADot = Dot()
        omegaPrime = (k/m)**0.5 * (1-(b**2)/(4*k*m))*0.5
        exponentialChange = Dot(color = RED)
        periodicChange = Dot(color = BLUE)
        time = 0
        d_t = 0.01
        ADot.shift(RIGHT*1.5, UP*1.2)

        combiningExplanation = TextMobject("A can be any random complex number").scale(0.7)
        combiningExplanation.move_to(functionUsed.get_center()+DOWN)
        self.play(Write(combiningExplanation), run_time = 2)
        self.wait()
        self.play(FadeOut(combiningExplanation))
        self.wait()

        plane = ScreenGrid()
        plane.shift(DOWN*0.5)
        self.play(ShowCreation(plane))
        self.play(Write(ADot))
        self.wait()

        reWritingA = TexMobject("A = |A|e^{i \\phi}").scale(0.5)
        magA = Line(start = DOWN*0.5, end = 1.5*RIGHT + 1.2*UP)
        magAName = TexMobject("|A|")
        angle = Arc(arc_center = DOWN*0.5, radius = 2.267, angle = np.arctan(1.7/1.5))
        angleName = TexMobject("\\phi")
        magAName.next_to(magA, LEFT, buff = 0)
        angleName.next_to(angle, RIGHT, buff = 0)
        reWritingFunction = TexMobject("x(t) = |A|e^{i \\phi}", " e^{\\frac{-b}{2m}t}", "e^{i\\omega\'t}")
        reWritingFunction[1].set_color(RED)
        reWritingFunction[2].set_color(BLUE)
        reWritingA.move_to(functionUsed.get_center()+ DOWN*0.5)
        reWritingFunction.to_edge(UP)

        self.play(Write(reWritingA), Write(magA), Write(angle),Write(magAName), Write(angleName), run_time = 4)
        self.wait()
        self.play(FadeOut(reWritingA), FadeOut(magA), FadeOut(angle), FadeOut(magAName), FadeOut(angleName))
        self.play(ReplacementTransform(functionUsed, reWritingFunction))
        self.wait()
        finalFunction = TexMobject("x(t)", " = |A|", "e^{\\frac{-b}{2m}t}", "e^{i(\\phi+\\omega\'t)}")
        finalFunction.to_edge(UP)
        finalFunction[0].set_color(GREEN)
        finalFunction[1].set_color(RED)
        finalFunction[2].set_color(RED)
        finalFunction[3].set_color(BLUE)
        self.play(ReplacementTransform(reWritingFunction, finalFunction))
        self.wait()

        displacement = Dot(color = GREEN)
        displacement.move_to(A*(np.exp(-b/(2*m) * time) * RIGHT * np.cos(omegaPrime*time + phi) + A*np.exp(-b/(2*m) * time) * np.sin(omegaPrime*time + phi)* UP + DOWN*0.5))
        exponentialChange.move_to(np.exp(-b/(2*m) * time) * RIGHT + DOWN*0.5)
        periodicChange.move_to(np.cos(omegaPrime*time + phi)*RIGHT + DOWN*0.5 + np.sin(omegaPrime*time + phi)* UP)
        self.play(Write(exponentialChange), Write(periodicChange))
        self.wait()
        time+=d_t

        for i in range(1000):
            newExponentialChange = Dot(color = RED)
            newExponentialChange.move_to(np.exp(-b/(2*m) * time) * RIGHT+ DOWN*0.5)
            newPeriodicChange = Dot(color = BLUE)
            newPeriodicChange.move_to(np.cos(omegaPrime*time + phi)*RIGHT + np.sin(omegaPrime*time + phi)* UP + DOWN*0.5)
            newDisplacement = Dot(color = GREEN)
            newDisplacement.move_to(A*(np.exp(-b/(2*m) * time) * RIGHT * np.cos(omegaPrime*time + phi) + A*np.exp(-b/(2*m) * time) * np.sin(omegaPrime*time + phi)* UP + DOWN*0.5))

            self.play(Transform(exponentialChange, newExponentialChange), Transform(periodicChange, newPeriodicChange), Transform(displacement, newDisplacement), run_time = 0.01)
            time+=d_t

        note = TextMobject("""The Green Dot shows how the box is supposed to move;
                                but it doesn't actually move like this!
        """).scale(0.7)
        self.play(ShowCreation(note))
        self.wait(2)



class ClearingThingsUp(Scene):
    def construct(self):
        finalFunction = TexMobject("x(t)", " = |A|", "e^{\\frac{-b}{2m}t}", "e^{i(\\phi+\\omega\'t)}")
        finalFunction.to_edge(UP)
        note = TextMobject("The box doesn't move upwards and downwards, along the imaginary line").scale(0.7)
        noteTwo = TextMobject("So there is no imaginary component - it can be ignored").scale(0.7)
        noteThree = TexMobject("e^{i(\\phi+\\omega\'t)} = cos(\\phi+\\omega\'t) + isin(\\phi+\\omega\'t)").scale(0.7)
        noteFour = TextMobject("Since there is no imaginary component, the imaginary sin component can be ignored!").scale(0.7)
        solution = TexMobject("x(t)=Ae^{\\frac{-b}{2m}t}cos(\\phi+\\omega\'t)")

        note.move_to(finalFunction.get_center()+DOWN*0.7)
        noteTwo.move_to(note.get_center()+DOWN*0.7)
        noteThree.move_to(noteTwo.get_center() + DOWN*0.7)
        noteFour.move_to(noteThree.get_center() + DOWN*0.7)

        self.play(FadeIn(finalFunction))
        self.wait(2)
        self.play(ShowCreation(note))
        self.wait(2)
        self.play(ShowCreation(noteTwo))
        self.wait(2)
        self.play(ShowCreation(noteThree))
        self.wait(2)
        self.play(ShowCreation(noteFour))
        self.wait(2)
        self.play(FadeOut(finalFunction), FadeOut(note), FadeOut(noteTwo), FadeOut(noteThree), FadeOut(noteFour))
        self.wait(2)
        self.play(Write(solution), run_time = 3)
        self.play(Indicate(solution))
        self.wait(2)

        A = TextMobject("A is the Amplitude of oscillation (where the box is let go from)").scale(0.5)
        b = TextMobject("b is the coefficient of dampening (how much friction/drag there is)").scale(0.5)
        m = TextMobject("m is the mass of the system").scale(0.5)
        phi = TexMobject("\\phi").scale(0.5)
        phiTwo = TextMobject(" is the phase of oscillation (what point on the wave the box is)").scale(0.5)
        omegaPrime = TexMobject("\\omega \' = \\sqrt{\\frac { k }{ m }}  \\sqrt { 1- \\frac { b^{ 2 } }{ 4km }  }").scale(0.5)
        omegaPrimeTwo = TextMobject(" can be divided by 2pi to find frequency of oscillation (how many waves are produced per second)").scale(0.5)


        A.move_to(solution.get_center()+DOWN*0.5)
        b.move_to(A.get_center()+DOWN*0.5)
        m.move_to(b.get_center()+DOWN*0.5)
        phi.move_to(m.get_center()+DOWN*0.5+LEFT*3)
        phiTwo.next_to(phi, RIGHT)
        omegaPrime.move_to(phi.get_center()+DOWN+LEFT*2)
        omegaPrimeTwo.next_to(omegaPrime)

        self.play(Write(A), Write(b), Write(m), Write(phi), Write(phiTwo), Write(omegaPrime), Write(omegaPrimeTwo))
        self.wait(8)
