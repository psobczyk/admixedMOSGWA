#include "QRuncher.hpp"
#include "TestSuite.hpp"
#include <math.h>

using namespace std;
using namespace linalg;
using namespace unitpp;

namespace test {

	/** Tests the class {@link QRuncher}. */
	class QRuncherTest : public TestSuite {

		/** Test with a sequence of actions. */
		void testSequence ();

		/** Test calculation of regression coefficients. */
		void testCoefficients ();

		/** Test the (default) copy constructor. */
		void testCopyConstructor ();

		/** Test the (default) assignment operator. */
		void testAssignment ();

		/** Tests {@link QRuncher::calculateSkipColumnRSS}. */
		void testSkipColumn ();

		/** Tests {@link QRuncher::calculateSkipColumnRSS}
		* in the case of skipping the last added linearly independent column.
		*/
		void testSkipLastLinearlyIndependentColumn ();

		public:

		QRuncherTest () : TestSuite( "QRuncherTest" ) {
			addTestMethod( "QRuncherTest::testSequence", this, &QRuncherTest::testSequence );
			addTestMethod( "QRuncherTest::testCoefficients", this, &QRuncherTest::testCoefficients );
			addTestMethod( "QRuncherTest::testCopyConstructor", this, &QRuncherTest::testCopyConstructor );
			addTestMethod( "QRuncherTest::testAssignment", this, &QRuncherTest::testAssignment );
			addTestMethod( "QRuncherTest::testSkipColumn", this, &QRuncherTest::testSkipColumn );
			addTestMethod( "QRuncherTest::testSkipLastLinearlyIndependentColumn", this, &QRuncherTest::testSkipLastLinearlyIndependentColumn );
		}
	} * qrUncherTest = new QRuncherTest();	// automatically freed by unit++

	void QRuncherTest::testSequence () {
		const double
			// col 1 norm 9, col 2 norm 7, col 3 = col 2 - col 1
			xData[] = {
				1, 2, 1,
				4, 3, -1,
				8, 6, -2
			},
			yData[] = { 1, 4, 8 };

		AutoMatrix x( 3, 3 );
		x.fill( xData );
		AutoVector y( 3 );
		y.fill( yData );

		QRuncher qr( y );
		assert_eq( "Zero variables", 81, qr.calculateRSS() );
		assert_eq( "Zero-dim regression coefficients", 0, qr.calculateCoefficients().countDimensions() );

		assert_eq(
			"No linear dep, full squares",
			9,
			fabs( qr.pushColumn( x.columnVector( 0 ) ) )
		);
		assert_close( "Exact match", 0, qr.calculateRSS() );
		assert_eq( "One-dim regression coefficients", 1, qr.calculateCoefficients().countDimensions() );
		assert_close( "Exact match coefficients", 1, qr.calculateCoefficients().get( 0 ) );

		AutoVector orthogonalPart( 3 );
		// col 1 * col 2 / norm^2 col1 = 62/81
		orthogonalPart.set( 0, 2 - 62./81. * 1 );
		orthogonalPart.set( 1, 3 - 62./81. * 4 );
		orthogonalPart.set( 2, 6 - 62./81. * 8 );
		assert_close(
			"No linear dep, reduced squares",
			sqrt( orthogonalPart.sumSquares() ),
			fabs( qr.pushColumn( x.columnVector( 1 ) ) ),
			1e-14
		);
		assert_close( "Still exact match", 0, qr.calculateRSS() );
		assert_eq( "Two-dim regression coefficients", 2, qr.calculateCoefficients().countDimensions() );
		assert_close( "Exact match coefficient 0", 1, qr.calculateCoefficients().get( 0 ) );
		assert_close( "Exact match coefficient 1", 0, qr.calculateCoefficients().get( 1 ) );

		assert_close(
			"Full linear dep, zero squares",
			0,
			fabs( qr.pushColumn( x.columnVector( 2 ) ) )
		);
		assert_close( "And still exact match", 0, qr.calculateRSS() );
		// Note: Regression coefficients are unspecified in linearly dependent case: therefore, no assertion

		qr.popColumn();
		assert_close( "Back to still exact match", 0, qr.calculateRSS() );
		assert_eq( "Back to two-dim regression coefficients", 2, qr.calculateCoefficients().countDimensions() );
		assert_close( "Back to exact match coefficient 0", 1, qr.calculateCoefficients().get( 0 ) );
		assert_close( "Back to exact match coefficient 1", 0, qr.calculateCoefficients().get( 1 ) );

		qr.popColumn();
		assert_close( "Back to exact match", 0, qr.calculateRSS() );
		assert_eq( "Back to one-dim regression coefficients", 1, qr.calculateCoefficients().countDimensions() );
		assert_close( "Back to exact match coefficients", 1, qr.calculateCoefficients().get( 0 ) );

		qr.popColumn();
		assert_eq( "Back to zero variables", 81, qr.calculateRSS() );
		assert_eq( "Back to zero-dim regression coefficients", 0, qr.calculateCoefficients().countDimensions() );
	}

	void QRuncherTest::testCoefficients () {
		const double
			xData[] = {
				3, 4,
				4, 3
			},
			yData[] = { 8, 8 };

		AutoMatrix x( 2, 2 );
		x.fill( xData );
		AutoVector y( 2 );
		y.fill( yData );

		QRuncher qr( y );
		assert_eq( "No variables yet", 0, qr.calculateCoefficients().countDimensions() );

		assert_eq( "No linear dep 0", 5, fabs( qr.pushColumn( x.columnVector( 0 ) ) ) );
		// Minimise || (8,8) - t * (3,4) ||^2 = (8-3t)^2+(8-4t)^2.
		// Derivative proportional 3(8-3t) + 4(8-4t) = 24-9t + 32-16t = 56-25t
		assert_close( "Approximate match coefficient 0", 56./25., qr.calculateCoefficients().get( 0 ) );

		assert_true( "No linear dep 1", 0.5 < fabs( qr.pushColumn( x.columnVector( 1 ) ) ) );
		assert_close( "Exact match", 0, qr.calculateRSS() );
		assert_close( "Exact match coefficients 0", 8./7., qr.calculateCoefficients().get( 0 ) );
		assert_close( "Exact match coefficients 1", 8./7., qr.calculateCoefficients().get( 1 ) );

		qr.popColumn();
		assert_close( "Back to approximate match coefficient 0", 56./25., qr.calculateCoefficients().get( 0 ) );
	}

	void QRuncherTest::testCopyConstructor () {

		// test vectors with norms: 15, 13, 25
		const double data[][3] = {
			{ 2, 5, 14 },
			{ 3, 4, 12 },
			{ 12, 15, 16 }
		};

		AutoVector u( 3 ), v( 3 ), w( 3 );
		u.fill( data[0] );
		v.fill( data[1] );
		w.fill( data[2] );

		QRuncher qr( u );
		QRuncher qrc0( qr );
		assert_eq( "Zero variables same RSS", qr.calculateRSS(), qrc0.calculateRSS() );
		assert_eq( "Zero variables same coefficients", qr.calculateCoefficients(), qrc0.calculateCoefficients() );

		const double lum0 = qr.pushColumn( v );
		assert_eq( "Copy still zero variables RSS", 225, qrc0.calculateRSS() );
		assert_eq( "Copy still zero variables coefficients", 0, qrc0.calculateCoefficients().countDimensions() );
		assert_eq( "Same linear dependence measure", lum0, qrc0.pushColumn( v ) );
		assert_eq( "Zero plus 1 variables same RSS", qr.calculateRSS(), qrc0.calculateRSS() );
		assert_eq( "Zero plus 1 variables same coefficients", qr.calculateCoefficients(), qrc0.calculateCoefficients() );

		QRuncher qrc1( qr );
		assert_eq( "One variable same RSS", qr.calculateRSS(), qrc1.calculateRSS() );
		assert_eq( "One variable same coefficients", qr.calculateCoefficients(), qrc1.calculateCoefficients() );

		const double lum1 = qr.pushColumn( w );
		assert_eq( "Copies still same one variable RSS", qrc0.calculateRSS(), qrc1.calculateRSS() );
		assert_eq( "Copies still same one variable coefficients", qrc0.calculateCoefficients(), qrc1.calculateCoefficients() );
		assert_eq( "Same linear dependence measure", lum1, qrc1.pushColumn( w ) );
		assert_eq( "Zero plus 1 variables same RSS", qr.calculateRSS(), qrc1.calculateRSS() );
		assert_eq( "Zero plus 1 variables same coefficients", qr.calculateCoefficients(), qrc1.calculateCoefficients() );

		QRuncher qrc2( qr );
		assert_eq( "Two variables same RSS", qr.calculateRSS(), qrc2.calculateRSS() );
		assert_eq( "Two variables same coefficients", qr.calculateCoefficients(), qrc2.calculateCoefficients() );

		qrc2.popColumn();
		assert_eq( "Same one var RSS after pop", qrc0.calculateRSS(), qrc2.calculateRSS() );
		assert_eq( "Same one var coefficients after pop", qrc0.calculateCoefficients(), qrc2.calculateCoefficients() );
		qrc2.popColumn();
		assert_eq( "Initial RSS after second pop", 225, qrc2.calculateRSS() );
	}

	void QRuncherTest::testAssignment () {

		// test vectors with norms: 17, 25, 29
		const double data[][2] = {
			{ 8, 15 },
			{ 7, 24 },
			{ 20, 21 }
		};

		AutoVector u( 2 ), v( 2 ), w( 2 );
		u.fill( data[0] );
		v.fill( data[1] );
		w.fill( data[2] );

		QRuncher qr( u );
		QRuncher qrc0( AutoVector( 0 ) );
		assert_eq( "Trivial RSS", 0, qrc0.calculateRSS() );
		assert_eq( "Trivial coefficients", 0, qrc0.calculateCoefficients().countDimensions() );

		qrc0 = qr;
		assert_eq( "Zero variables same RSS", qr.calculateRSS(), qrc0.calculateRSS() );
		assert_eq( "Zero variables same coefficients", qr.calculateCoefficients(), qrc0.calculateCoefficients() );

		const double lum0 = qr.pushColumn( v );
		assert_eq( "Copy still zero variables RSS", 289, qrc0.calculateRSS() );
		assert_eq( "Copy still zero variables coefficients", 0, qrc0.calculateCoefficients().countDimensions() );
		assert_eq( "Same linear dependence measure", lum0, qrc0.pushColumn( v ) );
		assert_eq( "Zero plus 1 variables same RSS", qr.calculateRSS(), qrc0.calculateRSS() );
		assert_eq( "Zero plus 1 variables same coefficients", qr.calculateCoefficients(), qrc0.calculateCoefficients() );

		QRuncher qrc1( v );
		qrc1 = qr;
		assert_eq( "One variable same RSS", qr.calculateRSS(), qrc1.calculateRSS() );
		assert_eq( "One variable same coefficients", qr.calculateCoefficients(), qrc1.calculateCoefficients() );

		const double lum1 = qr.pushColumn( w );
		assert_eq( "Copies still same one variable RSS", qrc0.calculateRSS(), qrc1.calculateRSS() );
		assert_eq( "Copies still same one variable coefficients", qrc0.calculateCoefficients(), qrc1.calculateCoefficients() );
		assert_eq( "Same linear dependence measure", lum1, qrc1.pushColumn( w ) );
		assert_eq( "Zero plus 1 variables same RSS", qr.calculateRSS(), qrc1.calculateRSS() );
		assert_eq( "Zero plus 1 variables same coefficients", qr.calculateCoefficients(), qrc1.calculateCoefficients() );

		QRuncher qrc2( u );
		qrc2 = qr;
		assert_eq( "Two variables same RSS", qr.calculateRSS(), qrc2.calculateRSS() );
		assert_eq( "Two variables same coefficients", qr.calculateCoefficients(), qrc2.calculateCoefficients() );

		qrc2.popColumn();
		assert_eq( "Same one var RSS after pop", qrc0.calculateRSS(), qrc2.calculateRSS() );
		assert_eq( "Same one var coefficients after pop", qrc0.calculateCoefficients(), qrc2.calculateCoefficients() );
		qrc2.popColumn();
		assert_eq( "Initial RSS after second pop", 289, qrc2.calculateRSS() );
	}

	void QRuncherTest::testSkipColumn () {
		const double
			// col norms: 11, 15, 18, 21, 18
			xData[] = {
				1,	2,	3,	-4,	-5,
				-2,	6,	5,	-5,	7,
				4,	8,	-11,	12,	-9,
				-10,	11,	-13,	16,	13
			},
			yData[] = { -2, -6, -8, -11 };

		AutoMatrix x( 4, 5 );
		x.fill( xData );
		AutoVector y( 4 );
		y.fill( yData );

		QRuncher qr( y );
		// Prepare comparison model RSS values
		qr.pushColumn( x.columnVector( 0 ) );
		qr.pushColumn( x.columnVector( 2 ) );
		double rss02 = qr.calculateRSS();
		qr.pushColumn( x.columnVector( 3 ) );
		double rss023 = qr.calculateRSS();
		qr.popColumn();
		assert_eq( "After pop is before push", rss02, qr.calculateRSS() );
		qr.popColumn();
		qr.popColumn();

		assert_eq( "Model 0: Zero variables, full RSS", 225, qr.calculateRSS() );

		assert_eq( "Model 1: No linear dep, full squares", 11, fabs( qr.pushColumn( x.columnVector( 0 ) ) ) );
		const double rss1 = qr.calculateRSS();
		assert_true( "Model 1: Partial match, significant remainder", 0.5 < rss1 );
		assert_eq( "Model 1: Back to zero variables", 225, qr.calculateSkipColumnRSS( 0 ) );
		assert_eq( "Model 1: RSS unchanged", rss1, qr.calculateRSS() );

		assert_true( "Model 2: No linear dep", 0.5 < fabs( qr.pushColumn( x.columnVector( 1 ) ) ) );
		const double rss2 = qr.calculateRSS();
		assert_close( "Model 2: Exact match", 0, rss2 );
		assert_close( "Model 2: Back to exact match", 0, qr.calculateSkipColumnRSS( 0 ) );
		assert_eq( "Model 2: Back to one variable, from y-stack", rss1, qr.calculateSkipColumnRSS( 1 ) );
		assert_eq( "Model2: RSS unchanged", rss2, qr.calculateRSS() );

		assert_true( "Model 3: No linear dep", 0.5 < fabs( qr.pushColumn( x.columnVector( 2 ) ) ) );
		const double rss3 = qr.calculateRSS();
		assert_close( "Model 3: Exact match", 0, rss3 );
		assert_close( "Model 3: Back to still exact match", 0, qr.calculateSkipColumnRSS( 0 ) );
		assert_close( "Model 3: Back to inexact match", rss02, qr.calculateSkipColumnRSS( 1 ), 1e-12 );
		assert_eq( "Model 3: Back to one variable, from y-stack", rss2, qr.calculateSkipColumnRSS( 2 ) );
		assert_eq( "Model 3: RSS unchanged", rss3, qr.calculateRSS() );

		assert_true( "Model 4: No linear dep", 0.5 < fabs( qr.pushColumn( x.columnVector( 3 ) ) ) );
		const double rss4 = qr.calculateRSS();
		assert_close( "Model 4: Exact match", 0, rss4 );
		assert_close( "Model 4: Back to still exact match", 0, qr.calculateSkipColumnRSS( 0 ) );
		assert_close( "Model 4: Back to inexact match", rss023, qr.calculateSkipColumnRSS( 1 ), 1e-12 );
		assert_close( "Model 4: Other back to still exact match", 0, qr.calculateSkipColumnRSS( 2 ) );
		assert_eq( "Model 4: Back to one variable, from y-stack", rss3, qr.calculateSkipColumnRSS( 3 ) );
		assert_eq( "Model 4: RSS unchanged", rss4, qr.calculateRSS() );

		assert_true( "Can pop", qr.popColumn() );

		// Note: In linearly dependent case, RSS may be understated
		assert_true( "Model 5: No linear dep", 0.5 < fabs( qr.pushColumn( x.columnVector( 4 ) ) ) );
		const double rss5 = qr.calculateRSS();
		assert_close( "Model 5: Exact match", 0, rss5 );
		assert_close( "Model 5: Back to still exact match", 0, qr.calculateSkipColumnRSS( 0 ) );
		assert_close( "Model 5: Again exact match", 0, qr.calculateSkipColumnRSS( 2 ) );
		assert_eq( "Model 5: Back to one variable, from y-stack", rss3, qr.calculateSkipColumnRSS( 3 ) );
		assert_eq( "Model 5: RSS unchanged", rss5, qr.calculateRSS() );
	}

	void QRuncherTest::testSkipLastLinearlyIndependentColumn () {
		const double
			xData[] = {
				1,	0,
				0,	1
			},
			yData[] = { 3, 4 };

		AutoMatrix x( 2, 2 );
		x.fill( xData );
		AutoVector y( 2 );
		y.fill( yData );

		QRuncher qr( y );
		assert_eq( "rss", 25.0, qr.calculateRSS() );
		assert_eq( "l.u.[0]", 1.0, fabs( qr.pushColumn( x.columnVector( 0 ) ) ) );
		assert_eq( "rss0", 16.0, qr.calculateRSS() );
		assert_eq( "rss0-0", 25.0, qr.calculateSkipColumnRSS( 0 ) );
		assert_eq( "l.u.[1]", 1.0, fabs( qr.pushColumn( x.columnVector( 1 ) ) ) );
		assert_eq( "rss01", 0.0, qr.calculateRSS() );
		assert_eq( "rss01-0", 9.0, qr.calculateSkipColumnRSS( 0 ) );
		assert_eq( "rss01-1", 16.0, qr.calculateSkipColumnRSS( 1 ) );
		assert_true( "pop 1", qr.popColumn() );
		assert_eq( "postpop rss0", 16.0, qr.calculateRSS() );
		assert_eq( "postpop rss0-0", 25.0, qr.calculateSkipColumnRSS( 0 ) );
		assert_true( "pop 0", qr.popColumn() );
		assert_eq( "postpop rss", 25.0, qr.calculateRSS() );
		assert_false( "pop too many", qr.popColumn() );
		assert_eq( "postpoptoomany rss", 25.0, qr.calculateRSS() );
	}

}
